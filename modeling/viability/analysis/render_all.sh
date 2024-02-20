mkdir -p nb_render

for f in *.py
do
jupytext --to notebook --execute $f --set-kernel python3
jupfn="${f%.*}".ipynb
cp -f $jupfn nb_render/$jupfn
rm $jupfn
done