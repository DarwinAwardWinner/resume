#!/usr/bin/env doit -f

import regex
from fnmatch import fnmatch
import os.path
from subprocess import check_output
from doit import create_after
from distutils.spawn import find_executable
try:
    from os import scandir, walk
except ImportError:
    from scandir import scandir, walk

DOIT_CONFIG = {
    'default_tasks': ['publish_to_mneme'],
}

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

LYXPATH = find_executable("lyx") or "/Applications/LyX.app/Contents/MacOS/lyx"

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
            'file_dep': [lyxfile],
            'targets': [pdffile],
            'verbosity': 0,
        }

def task_readme2index():
    yield {
        'name': None,
        'doc': "Convert README.mkdn file to index.html."
    }
    for mkdnfile in glob_recursive("README.mkdn", top="examples"):
        htmlfile = os.path.join(os.path.dirname(mkdnfile), "index.html")
        yield {
            'name': mkdnfile,
            'actions': [["pandoc", "-t", "html", "-o", htmlfile, mkdnfile]],
            'file_dep': [mkdnfile],
            'targets': [htmlfile],
        }

@create_after(executed='readme2index')
def task_publish_to_mneme():
    rsync_common_args = [
        "rsync", "-vrL", "--size-only", "--delete",
        "--exclude", ".DS_Store", "--delete-excluded",
    ]
    rsync_srcs = ["ryan_thompson_resume.pdf", "examples"]
    rsync_dest = "mneme:public_html/resume/"
    # Compute file dependencies by running "rsync --list-only" and
    # parsing the output.
    rsync_list_cmd = rsync_common_args + rsync_srcs + [".", "--list-only"]
    rsync_out = check_output(rsync_list_cmd).splitlines()
    file_deps = []
    for line in rsync_out:
        s = regex.search("^-(?:\S+\s+){4}(.*)", line)
        if s is not None:
            file_deps.append(s.group(1))
    rsync_xfer_cmd = rsync_common_args + rsync_srcs + [rsync_dest]
    return {
        'actions': [rsync_xfer_cmd],
        'file_dep': file_deps,
        'doc': "Sync resume and supporting files to mneme.",
        'verbosity': 2,
    }
