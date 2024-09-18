COL_NAME=$1
INFILE=$2

awk -F'\t' -v col_name="$COL_NAME" '
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i == col_name) {
                col=i;
                break;
            }
        }
    }
    col {print $col}
' $INFILE
