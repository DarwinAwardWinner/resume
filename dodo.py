#!/usr/bin/env doit -f
# -*- coding: utf-8; -*-

from collections import Iterable, Mapping  # in Python 3 use from collections.abc
from distutils.spawn import find_executable
from fnmatch import fnmatch
from kwonly import kwonly
from subprocess import check_output
from types import StringTypes
import os.path
import regex

try:
    from os import scandir, walk
except ImportError:
    from scandir import scandir, walk

DOIT_CONFIG = {
    'default_tasks': ['publish_to_mneme'],
}

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
        if isinstance(arg, StringTypes):
            # String
            result.append(arg)
        elif isinstance(arg, Mapping):
            # Dict-like
            result.append(arg)
        elif isinstance(arg, Iterable):
            try:
                # Duck-typing test for list-ness (a stricter condition
                # that just "iterable")
                for i in xrange(len(arg)):
                    result.append(arg[i])
            except TypeError:
                # Not list-like
                result.append(arg)
        else:
            result.append(arg)
    return result

def find_mac_app(name):
    try:
        return check_output(
            ["mdfind",
             "kMDItemDisplayName==%s&&kMDItemKind==Application" % (name,) ]).strip()
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

@kwonly(0)
def rsync_list_files(extra_rsync_args=(), include_dirs=False, *paths):
    """Iterate over the files in path that rsync would copy.

By default, only files are listed, not directories, since doit doesn't
like dependencies on directories because it can't hash them.

This uses "rsync --list-only" to make rsync directly indicate which
files it would copy, so any exclusion/inclusion rules are taken into
account.

    """
    rsync_list_cmd = [ 'rsync', '-r', "--list-only" ] + unnest(extra_rsync_args) + unnest(paths) + [ "." ]
    rsync_out = check_output(rsync_list_cmd).splitlines()
    for line in rsync_out:
        s = regex.search("^(-|d)(?:\S+\s+){4}(.*)", line)
        if s is not None:
            if include_dirs or s.group(1) == '-':
                yield s.group(2)

rsync_common_args = ["-rL", "--size-only", "--delete", "--exclude", ".DS_Store", "--delete-excluded",]

def task_lyx2pdf():
    yield {
        'name': None,
        'doc': "Convert LyX file to PDF."
    }
    for lyxfile in glob_recursive("*.lyx"):
        pdffile = lyxfile[:-3] + "pdf"
        lyx_cmd = [LYXPATH, "--export-to", "pdf4" , pdffile, lyxfile]
        yield {
            'name': lyxfile,
            'actions': [lyx_cmd],
            # Assume every bib file is a dependency of any LaTeX
            # operation
            'file_dep': [lyxfile] + list(glob_recursive('*.bib')),
            'targets': [pdffile],
            'verbosity': 0,
            'clean': True,
        }

def task_readme2index():
    yield {
        'name': None,
        'doc': "Convert README.mkdn files to index.html."
    }
    for mkdnfile in glob_recursive("README.mkdn", top="examples"):
        htmlfile = os.path.join(os.path.dirname(mkdnfile), "index.html")
        yield {
            'name': mkdnfile,
            'actions': [["pandoc", "-t", "html", "-o", htmlfile, mkdnfile]],
            'file_dep': [mkdnfile],
            'targets': [htmlfile],
            'clean': True,
        }

def task_publish_to_mneme():
    yield {
        'name': None,
        'doc': "Sync resume and supporting files to mneme."
    }
    # Resume PDF file
    rsync_src = "ryan_thompson_resume.pdf"
    rsync_dest = "mneme:public_html/resume/"
    file_deps = [rsync_src]
    rsync_xfer_cmd = ["rsync"] + rsync_common_args + [ rsync_src, rsync_dest ]
    yield {
        'name': rsync_src,
        'actions': [rsync_xfer_cmd],
        'file_dep': file_deps,
        'task_dep': [ 'lyx2pdf' ],
        'doc': "rsync resume PDF file to mneme.",
        'verbosity': 2,
    }
    # Examples directory
    rsync_src = "examples"
    rsync_dest = "mneme:public_html/resume/"
    file_deps = list(rsync_list_files(rsync_src, extra_rsync_args=rsync_common_args))
    # Ensure the generated html files are in file_deps
    for f in file_deps:
        if os.path.basename(f) == "README.mkdn":
            htmlfile = os.path.join(os.path.dirname(f), "index.html")
            if htmlfile not in file_deps:
                file_deps.append(htmlfile)
    file_deps = sorted(file_deps)
    rsync_xfer_cmd = ["rsync"] + rsync_common_args + [ rsync_src, rsync_dest ]
    yield {
        'name': rsync_src,
        'actions': [rsync_xfer_cmd],
        'file_dep': file_deps,
        'task_dep': [ 'readme2index' ],
        'doc': "rsync examples directory to mneme.",
        'verbosity': 2,
    }

# Dummy target if you just want to build everything but not publish
def task_build():
    return {
        'doc': 'Build resume and supporting files.',
        'task_dep': [ 'lyx2pdf', 'readme2index' ],
        'actions': [],
    }
