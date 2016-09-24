# -*- coding: utf-8; -*-

import locale
import os.path
import regex

from collections import Iterable, Mapping  # in Python 3 use from collections.abc
from distutils.spawn import find_executable
from fnmatch import fnmatch
from subprocess import check_output, check_call
from tempfile import NamedTemporaryFile

try:
    from os import scandir, walk
except ImportError:
    from scandir import scandir, walk

def unnest(*args):
    """Un-nest list- and tuple-like elements in arguments.

"List-like" means anything with a len() and whose elments can be
accessed with numeric indexing, except for string-like elements. It
must also be an instance of the collections.Iterable abstract class.
Dict-like elements and iterators/generators are not affected.

This function always returns a list, even if it is passed a single
scalar argument.

    """
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
    """Shortcut for check.output + str.decode"""
    return check_output(*args, **kwargs).decode(encoding)

def find_mac_app(name):
    try:
        return check_output_decode(
            ["mdfind",
             "kMDItemDisplayName=={name}&&kMDItemKind==Application".format(name=name)]).split("\n")[0]
    except Exception:
        return None

def glob_recursive(pattern, top=".", include_hidden=False, *args, **kwargs):
    """Combination of glob.glob and os.walk.

Reutrns the relative path to every file or directory matching the
pattern anywhere in the specified directory hierarchy. Defaults to the
current working directory. Any additional arguments are passed to
os.walk."""
    for (path, dirs, files) in walk(top, *args, **kwargs):
        for f in dirs + files:
            if include_hidden or f.startswith("."):
                continue
            if fnmatch(f, pattern):
                yield os.path.normpath(os.path.join(path, f))

LYXPATH = find_executable("lyx") or \
    os.path.join(find_mac_app("LyX"), "Contents/MacOS/lyx") or \
    '/bin/false'

def rsync_list_files(*paths, extra_rsync_args=(), include_dirs=False):
    """Iterate over the files in path that rsync would copy.

By default, only files are listed, not directories, since doit doesn't
like dependencies on directories because it can't hash them.

This uses "rsync --list-only" to make rsync directly indicate which
files it would copy, so any exclusion/inclusion rules are taken into
account.

    """
    rsync_list_cmd = [ 'rsync', '-r', "--list-only" ] + unnest(extra_rsync_args) + unnest(paths) + [ "." ]
    rsync_out = check_output_decode(rsync_list_cmd).splitlines()
    for line in rsync_out:
        s = regex.search("^(-|d)(?:\S+\s+){4}(.*)", line)
        if s is not None:
            if include_dirs or s.group(1) == '-':
                yield s.group(2)

def lyx_image_deps(wildcards):
    lyxfile = wildcards.filename + ".lyx"


def lyx_bib_deps(wildcards):
    # Cheat: Assume every bib file is a dependency of any LaTeX
    # operation
    return list(glob_recursive('*.bib'))

readme_files = list(glob_recursive("README.mkdn", top="examples"))
index_files = [ os.path.join(os.path.dirname(f), "index.html") for f in readme_files ]

rsync_common_args = ["-rL", "--size-only", "--delete", "--exclude", ".DS_Store", "--delete-excluded",]

all_example_files = set(rsync_list_files('examples', extra_rsync_args=rsync_common_args))
all_example_files = all_example_files.union(index_files)

rsync_dest = "mneme:public_html/resume/"

rule build_all:
    input: "ryan_thompson_resume.pdf", index_files

rule create_resume_pdf:
    input: lyxfile="ryan_thompson_resume.lyx", bibfile="citations.bib", headshot="headshot-crop.jpg"
    output: pdf="ryan_thompson_resume.pdf"
    shell: '{LYXPATH:q} --export-to pdf4 {output.pdf:q} {input.lyxfile:q}'

rule create_resume_html:
    input: lyxfile="ryan_thompson_resume.lyx", bibfile="citations.bib", headshot="headshot-crop.jpg"
    output: html="ryan_thompson_resume.html"
    run:
        with NamedTemporaryFile() as tempf:
            shell('{LYXPATH:q} --export-to xhtml {tempf.name:q} {input.lyxfile:q}')
            shell('''cat {tempf.name:q} | perl -lape 's[<span class="flex_cv_image">(.*?)</span>][<span class="flex_cv_image"><img src="$1" width="100"></span>]g' > {output.html:q}''')

rule link_resume_to_index_html:
    input: 'ryan_thompson_resume.html'
    output: 'index.html'
    shell: 'ln -s {input:q} {output:q}'

rule readme_to_index_html:
    input: "{dirname}/README.mkdn"
    output: "{dirname}/index.html"
    shell: 'pandoc -t html -o {output[0]:q} {input[0]:q}'

rule publish:
    input: roots=('ryan_thompson_resume.pdf', 'ryan_thompson_resume.html', 'index.html', 'examples', 'headshot-crop.jpg'),
           others=all_example_files
    shell: '''
    rsync --info=progress2 {rsync_common_args:q} \
      {input.roots:q} {rsync_dest:q}
    '''
