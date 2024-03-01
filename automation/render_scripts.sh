: '
script to render all python scripts in the current directory to notebooks
'

# use -f to force re-rendering of all notebooks
force=0
while getopts "f" opt; do
  case ${opt} in
    f)
      force=1
      ;;
    \?)
      echo "Invalid option: $OPTARG" 1>&2
      ;;
  esac
done
shift $((OPTIND -1))

mkdir -p nb_render

for f in *.py
do
  jupfn="nb_render/${f%.*}.ipynb"
  if [[ ! -f "$jupfn" || $force -eq 1 ]]; then
    echo "---Rendering $f to $jupfn---"
    jupytext --to notebook --execute $f --set-kernel python3
    mv "${f%.*}.ipynb" "$jupfn"
  else
    echo "---Skipping $f as notebook already exists---"
  fi
done