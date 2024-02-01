Cannot get markdown snippets to render with markdown package

Appears to be because some render markdown.in file is being put in the output directory, which can then not be read in a subsequent step. 

Giving up on files in output directory. going to try workflow of cleaning files 

https://tex.stackexchange.com/a/87846

Setup a recie that runs latexmk with lualatex then runs with the -c option to clean. Doesn't seem to delete all extra files, but just gitignoring those. 

useful
https://www.joelotz.com/blog/2022/enabling-different-latex-compilers-in-vs-code.html

https://applegamer22.github.io/posts/tex/vscode/



## citations

using zotero export with better bibtex. use 'export' in zotero to .bib file. 

https://retorque.re/zotero-better-bibtex/citing/cayw/index.html

Need to do /cite{Citation Key}

There is a vscode extension for not having to type the citation key, but appears to not work on WSL

https://github.com/mblode/vscode-zotero

https://github.com/mblode/vscode-zotero/issues/14