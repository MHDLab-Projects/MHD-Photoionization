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


## word conversion 

### Word to latex

https://tanghaoyu258.github.io/posts/2021/04/blog-post-1/

`pandoc --extract-media=. --wrap=none -s input.docx -t latex –o output.tex`


### Latex to word

https://tex.stackexchange.com/questions/111886/how-to-convert-a-scientific-manuscript-from-latex-to-word-using-pandoc

https://latex2rtf.sourceforge.net/

`latex2rtf main.tex`

pandoc: `pandoc main.tex -o test.docx`

Just open pdf in word. 

## Interdocument references

use the xr package for references

The .aux files need to be compiled and kept for the references to work. Need to add \externaldocument for the specific document that the label is in. 

when using \ref, the reference is a hyperlink to the figure with the same number in the document...just going to use non hyperlink references with the SI using \ref*