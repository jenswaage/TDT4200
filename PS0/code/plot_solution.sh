help()
{
    echo
    echo "Plot 1D shallow water equations" 
    echo
    echo "Syntax"
    echo "---------------------------------------------"
    echo "./plot_solution.sh [-n|h]                    "
    echo 
    echo "Option    Description     Arguments   Default"
    echo "---------------------------------------------"
    echo "n         Grid size       Optional    1024   "
    echo "h         Help            None               "
    echo
    echo "Example"
    echo "---------------------------------------------"
    echo "./plot_solution.sh -n 1024                   "
    echo 
}

#-----------------------------------------------------------------
set -e

N=1024

while getopts ":n:h" opt; do
    case $opt in
        n)
            N=$OPTARG;;
        h)
            help
            exit;;
        \?)
            echo "Invalid option"
            help
            exit;;
    esac
done

#-----------------------------------------------------------------
SIZE=`echo $N | bc`

for FILE in `ls data/*.bin`
do
OUTFILE=`echo $FILE | sed s/^data/plots/ |  sed s/\.bin/.png/`
cat <<END_OF_SCRIPT | gnuplot -
set term png
set output "$OUTFILE"
set xrange [0:$SIZE]
set yrange [0:1]
plot "$FILE" binary format='%double' using 0:1
END_OF_SCRIPT
echo $OUTFILE
done
