# Run all python scripts and echo their name to the terminal

for i in $(ls *.py); do
    echo $i
    python $i
done