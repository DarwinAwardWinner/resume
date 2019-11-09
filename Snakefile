# -*- coding: utf-8; -*-

import locale
import os.path
import regex
import urllib.parse
import bibtexparser

from collections.abc import Iterable, Mapping
from distutils.spawn import find_executable
from fnmatch import fnmatch
from subprocess import check_output, check_call
from tempfile import NamedTemporaryFile
from bibtexparser.bibdatabase import BibDatabase

try:
    from os import scandir, walk
except ImportError:
    from scandir import scandir, walk

def unnest(*args):
    '''Un-nest list- and tuple-like elements in arguments.

"List-like" means anything with a len() and whose elments can be
accessed with numeric indexing, except for string-like elements. It
must also be an instance of the collections.Iterable abstract class.
Dict-like elements and iterators/generators are not affected.

This function always returns a list, even if it is passed a single
scalar argument.

    '''
    result = []
    for arg in args:
        if isinstance(arg, str):
            # String
            result.append(arg)
        elif isinstance(arg, Mapping):
            # Dict-like
            result.append(arg)
        elif isinstance(arg, Iterable):
            try:
                # Duck-typing test for list-ness (a stricter condition
                # than just "iterable")
                for i in range(len(arg)):
                    result.append(arg[i])
            except TypeError:
                # Iterable but not list-like
                result.append(arg)
        else:
            # Not iterable
            result.append(arg)
    return result

def check_output_decode(*args, encoding=locale.getpreferredencoding(), **kwargs):
    '''Shortcut for check.output + str.decode'''
    return check_output(*args, **kwargs).decode(encoding)

def find_mac_app(name):
    try:
        result = check_output_decode(
            ['mdfind',
             'kMDItemDisplayName=={name}.app&&kMDItemKind==Application'.format(name=name)]).split('\n')[0]
        if not result:
            raise Exception("No result found")
        return result
    except Exception:
        return None

def glob_recursive(pattern, top='.', include_hidden=False, *args, **kwargs):
    '''Combination of glob.glob and os.walk.

Reutrns the relative path to every file or directory matching the
pattern anywhere in the specified directory hierarchy. Defaults to the
current working directory. Any additional arguments are passed to
os.walk.'''
    for (path, dirs, files) in walk(top, *args, **kwargs):
        for f in dirs + files:
            if include_hidden or f.startswith('.'):
                continue
            if fnmatch(f, pattern):
                yield os.path.normpath(os.path.join(path, f))

LYXPATH = find_executable('lyx') or \
    os.path.join(find_mac_app('LyX'), 'Contents/MacOS/lyx') or \
    '/bin/false'

def rsync_list_files(*paths, extra_rsync_args=(), include_dirs=False):
    '''Iterate over the files in path that rsync would copy.

By default, only files are listed, not directories, since doit doesn't
like dependencies on directories because it can't hash them.

This uses "rsync --list-only" to make rsync directly indicate which
files it would copy, so any exclusion/inclusion rules are taken into
account.

    '''
    rsync_list_cmd = [ 'rsync', '-r', '--list-only' ] + unnest(extra_rsync_args) + unnest(paths) + [ '.' ]
    rsync_out = check_output_decode(rsync_list_cmd).splitlines()
    for line in rsync_out:
        s = regex.search('^(-|d)(?:\S+\s+){4}(.*)', line)
        if s is not None:
            if include_dirs or s.group(1) == '-':
                yield s.group(2)

def lyx_bib_deps(lyxfile):
    '''Return an iterator over all bib files referenced by a Lyx file.

    This will only return the names of existing files, so it will be
    unreliable in the case of an auto-generated bib file.

    '''
    try:
        with open(lyxfile) as f:
            lyx_text = f.read()
        bib_names = []
        for m in regex.finditer('bibfiles "(.*?)"', lyx_text):
            bib_names.extend(m.group(1).split(','))
        # Unfortunately LyX doesn't indicate which bib names refer to
        # files in the current directory and which don't. Currently that's
        # not a problem for me since all my refs are in bib files in the
        # current directory.
        for bn in bib_names:
            bib_path = bn + '.bib'
            yield bib_path
    except FileNotFoundError:
        pass

def lyx_hrefs(lyxfile):
    '''Return an iterator over hrefs in a LyX file.'''
    pattern = '''
    (?xsm)
    ^ LatexCommand \\s+ href \\s* \\n
    (?: name \\b [^\\n]+ \\n )?
    target \\s+ "(.*?)" $
    '''
    with open(lyxfile) as f:
        return (urllib.parse.unquote(m.group(1)) for m in
                re.finditer(pattern, f.read()))

examples_base_url = 'https://darwinawardwinner.github.io/resume/examples/'
examples_dir = 'examples'

def resume_example_deps(lyxfile):
    '''Iterate over all referenced example files in a LyX file.'''
    for href in lyx_hrefs(lyxfile):
        if href.startswith(examples_base_url) and not href.endswith('/'):
            expath = href[len(examples_base_url):]
            yield os.path.join(examples_dir, expath)

readme_files = list(glob_recursive('README.mkdn', top='examples'))
index_files = [ os.path.join(os.path.dirname(f), 'index.html') for f in readme_files ]

rsync_common_args = ['-rL', '--size-only', '--delete', '--exclude', '.DS_Store', '--delete-excluded',]

all_example_files = set(rsync_list_files('examples', extra_rsync_args=rsync_common_args))
r_html_files = [ f + '.html' for f in all_example_files if f.endswith('.R') ]
all_example_files = all_example_files.union(index_files)
all_example_files = all_example_files.union(r_html_files)

rule build_all:
    input: 'ryan_thompson_resume.pdf', 'ryan_thompson_resume.html', 'index.html', index_files, r_html_files

rule create_resume_pdf:
    input: lyxfile='ryan_thompson_resume.lyx',
           bibfiles=list(lyx_bib_deps('ryan_thompson_resume.lyx')),
           example_files=list(resume_example_deps('ryan_thompson_resume.lyx')),
           headshot='headshot-crop.png',
    output: pdf='ryan_thompson_resume.pdf'
    shell: '{LYXPATH:q} --export-to pdf4 {output.pdf:q} {input.lyxfile:q}'

rule create_resume_html:
    input: lyxfile='ryan_thompson_resume.lyx',
           bibfiles=list(lyx_bib_deps('ryan_thompson_resume.lyx')),
           example_files=list(resume_example_deps('ryan_thompson_resume.lyx')),
           headshot='headshot-crop.png',
    output: html='ryan_thompson_resume.html'
    run:
        with NamedTemporaryFile() as tempf:
            shell('{LYXPATH:q} --export-to xhtml {tempf.name:q} {input.lyxfile:q}')
            shell('''cat {tempf.name:q} | perl -lape 's[<span class="flex_cv_image">(.*?)</span>][<span class="flex_cv_image"><img src="$1" width="100"></span>]g' > {output.html:q}''')

rule create_cv_pdf:
    input: lyxfile='ryan_thompson_cv.lyx',
           bibfiles=list(lyx_bib_deps('ryan_thompson_cv.lyx')),
           example_files=list(resume_example_deps('ryan_thompson_cv.lyx')),
           headshot='headshot-crop.png',
    output: pdf='ryan_thompson_cv.pdf'
    shell: '{LYXPATH:q} --export-to pdf4 {output.pdf:q} {input.lyxfile:q}'

rule create_cv_html:
    input: lyxfile='ryan_thompson_cv.lyx',
           bibfiles=list(lyx_bib_deps('ryan_thompson_cv.lyx')),
           example_files=list(resume_example_deps('ryan_thompson_cv.lyx')),
           headshot='headshot-crop.png',
    output: html='ryan_thompson_cv.html'
    run:
        with NamedTemporaryFile() as tempf:
            shell('{LYXPATH:q} --export-to xhtml {tempf.name:q} {input.lyxfile:q}')
            shell('''cat {tempf.name:q} | perl -lape 's[<span class="flex_cv_image">(.*?)</span>][<span class="flex_cv_image"><img src="$1" width="100"></span>]g' > {output.html:q}''')

rule link_resume_to_index_html:
    input: 'ryan_thompson_resume.html'
    output: 'index.html'
    shell: 'ln -s {input:q} {output:q}'

rule examples_readme_to_index_html:
    input: '{dirname}README.mkdn'
    output: '{dirname,examples(/.*)?/}index.html'
    shell: 'pandoc -t html -o {output[0]:q} {input[0]:q}'

rule R_to_html:
    input: '{dirname}/{basename,[^/]+}.R'
    output: '{dirname}/{basename}.R.html'
    shell: 'pygmentize -f html -O full -l R -o {output:q} {input:q}'

rule process_bib:
    '''Preprocess bib file for LaTeX.

For entries with a DOI, all URLs are stripped, since the DOI already
provides a clickable link. For entries with no DOI, all but one URL is
discarded, since LyX can't handle entries with multiple URLs. The
shortest URL is kept.'''
    input: '{basename}.bib'
    output: '{basename,.*(?<!-PROCESSED)}-PROCESSED.bib'
    run:
        with open(input[0]) as infile:
            bib_db = bibtexparser.load(infile)
        entries = bib_db.entries
        for entry in entries:
            # Keep DOI or exactly one URL
            if 'doi' in entry:
                try:
                    del entry['url']
                except KeyError:
                    pass
            else:
                try:
                    entry_urls = regex.split('\\s+', entry['url'])
                    shortest_url = min(entry_urls, key=len)
                    # Need to fix e.g. 'http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=55329{\\&}tool=pmcentrez{\\&}rendertype=abstract'
                    shortest_url = re.sub('\\{\\\\?(.)\\}', '\\1', shortest_url)
                    entry['url'] = shortest_url
                except KeyError:
                    pass
            # Boldface my name
            authors = regex.split("\\s+and\\s+",entry['author'])
            for i in range(len(authors)):
                m = regex.search('^Thompson,\\s+(R.*)$', authors[i])
                if m:
                    authors[i] = f"\\textbf{{{m.group(1)} Thompson}}"
            entry['author'] = ' and '.join(authors)
        new_db = BibDatabase()
        new_db.entries = entries
        with open(output[0], 'w') as outfile:
            bibtexparser.dump(new_db, outfile)
