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

# Create or clear the error log file
echo "" > nb_render/error_log.txt

for f in *.py
do
  jupfn="nb_render/${f%.*}.ipynb"
  if [[ ! -f "$jupfn" || $force -eq 1 ]]; then
    echo "---Rendering $f to $jupfn---"
    echo -e "\n   - - -   Rendering $f to $jupfn   - - -\n" >> nb_render/error_log.txt
    # Redirect stderr to a temporary file
    jupytext --to notebook --execute $f --set-kernel python3 2> nb_render/temp_err.txt
    exit_status=$?
    if [ $exit_status -ne 0 ]; then
      echo "Error: $f failed to execute." >> nb_render/error_log.txt
      # Append the error message to the error log file
      # Append the error message to the error log file
      # Use sed to remove ANSI escape codes
      sed 's/\x1b\[[0-9;]*m//g' nb_render/temp_err.txt >> nb_render/error_log.txt
    else
      mv "${f%.*}.ipynb" "$jupfn"
    fi
  else
    echo "---Skipping $f as notebook already exists---"
  fi
done
