This is my résumé, built using Rob Oakes' excellent
[xetexCV class](https://code.google.com/archive/p/latex-professional/)
for [LaTeX](https://www.latex-project.org/) using 
[LyX](https://www.lyx.org/), along with a series of examples of
my work that are linked to from the résumé itself. I use
the [Snakemake](https://snakemake.readthedocs.io/en/stable/) build
tool to automate generating the PDF of my résumé, preparing all the
supporting files, and deploying them all to a web server for online
viewing. You can see the result here:
https://darwinawardwinner.github.io/resume/ryan_thompson_resume.pdf

Lastly, there is a manually-prepared plain text version here:
https://raw.githubusercontent.com/DarwinAwardWinner/resume/master/ryan_thompson_resume.txt

If you want to use this as a template for your own résumé, you'll need
to install a few fonts (or else switch the style to fonts that you
prefer):

1. Minion Pro: https://typekit.com/fonts/minion-pro
2. Fontin: https://www.exljbris.com/fontin.html
3. Fontin Sans: https://www.exljbris.com/fontinsans.html

Additionally, my résumé uses slightly modified versions of the xetexCV
document class ([`xetexCV.cls`](./xetexCV.cls)) and LyX layout
([`xetexCV.layout`](./xetexCV.layout)), which are included in this
repository. You can install these files to the appropriate locations
for
[LaTeX](http://blog.pengyifan.com/where-to-place-you-own-sty-or-cls-files-to-make-them-available-to-all-my-tex-files/)
and [LyX](https://wiki.lyx.org/Layouts/Layouts), but they should also
work if you just place them in the same directory as the LyX file.

## License

All files in this repository are licensed under the [Creative Commons
Attribution-ShareAlike 4.0
license](https://creativecommons.org/licenses/by-sa/4.0/), with the
exception of `xetexCV.cls` and `xetexCV.layout`. These two files are
based on the originals created by Rob Oakes, which are licensed under
the [LGPL license](https://www.gnu.org/licenses/lgpl-3.0.en.html).
Hence, my modified versions are also licensed under the same terms.
(Rob Oakes' [original 2009 mailing list
post](https://tug.org/pipermail/xetex/2009-December/015046.html)
doesn't specify a version of the LGPL, but the latest version at that
time was version 3, so I'm assuming that's what was meant.)
