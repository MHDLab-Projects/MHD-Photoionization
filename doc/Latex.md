Cannot get markdown snippets to render with markdown package

Appears to be because some render markdown.in file is being put in the output directory, which can then not be read in a subsequent step. 

Giving up on files in output directory. going to try workflow of cleaning files 

https://tex.stackexchange.com/a/87846

Setup a recie that runs latexmk with lualatex then runs with the -c option to clean. Doesn't seem to delete all extra files, but just gitignoring those. 

useful
https://www.joelotz.com/blog/2022/enabling-different-latex-compilers-in-vs-code.html

https://applegamer22.github.io/posts/tex/vscode/
