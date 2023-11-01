for f in analysis*.py
do
jupytext --to notebook --execute $f --set-kernel python3
jupfn="${f%.*}".ipynb
cp -f $jupfn nb_render/$jupfn
rm $jupfn
done