#!/bin/bash

# print usage
usage() {
    echo "Usage: $0 <binary> [-f <flag_file> | -np <NP> -p <NR> -q <NC> -NB <QC> [-MB <QR>] [--debug] [--use-scr-need-checkpoint]]"
    exit 1
}

# Verify if a string is an integer
is_integer() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

# Process the flags from a file
process_flags_from_file() {
    local unknown_flags=()

    while read -r line || [[ -n "$line" ]]; do
        # Verifica se a linha é vazia ou um comentário
        if [[ -z "$line" || "$line" =~ ^\s*# ]]; then
            continue
        fi

        # Extrai a flag e seu valor
        flag=$(echo "$line" | awk '{print $1}')
        value=$(echo "$line" | awk '{print $2}')

        case "$flag" in
            -np)
                np=$value
                ;;
            -p)
                nr=$value
                if ! is_integer "$nr"; then
                    echo "The value of -p must be an integer."
                    usage
                fi
                ;;
            -q)
                nc=$value
                if ! is_integer "$nc"; then
                    echo "The value of -q must be an integer."
                    usage
                fi
                ;;
            -NB)
                qc=$value
                if ! is_integer "$qc"; then
                    echo "The value of -NB must be an integer."
                    usage
                fi
                ;;
            -MB)
                mb=$value
                if ! is_integer "$mb"; then
                    echo "The value of -MB must be an integer."
                    usage
                fi
                ;;
            --debug)
                debug=true
                ;;
            --use-scr-need-checkpoint)
                use_scr_need_checkpoint=true
                ;;
            *)
                unknown_flags+=("$flag")
                ;;
        esac
    done < "$1"

    # Verify whether there are unknown flags and print them
    if [ ${#unknown_flags[@]} -gt 0 ]; then
        for flag in "${unknown_flags[@]}"; do
            echo "Unknown flag: $flag"
        done
        usage
    fi
}

# Define the variables
binary=$1
shift # move to the next argument

# Verify whether the binary was provided and is executable
if [ ! -x "./$binary" ]; then
    echo "The binary '$binary' was not found or is not executable."
    exit 1
fi

file=""
np=""
nr=""
nc=""
qc=""
qr=""
debug=false
use_scr_need_checkpoint=false
mb=""

# Verify whether the arguments were provided from a file or from the command line
if [ "$1" == "-f" ]; then
    if [ "$#" -lt 2 ]; then
        usage
    fi
    file="$2"
    if [ ! -f "$file" ]; then
        echo "Arquivo de flags não encontrado: $file"
        usage
    fi
    process_flags_from_file "$file"
else
    # Extract the arguments
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -np)
                np=$2
                shift
                shift
                ;;
            -p)
                nr=$2
                if ! is_integer "$nr"; then
                    echo "The value of -p must be an integer."
                    usage
                fi
                shift
                shift
                ;;
            -q)
                nc=$2
                if ! is_integer "$nc"; then
                    echo "The value of -q must be an integer."
                    usage
                fi
                shift
                shift
                ;;
            -NB)
                qc=$2
                if ! is_integer "$qc"; then
                    echo "The value of -NB must be an integer."
                    usage
                fi
                shift
                shift
                ;;
            -MB)
                mb=$2
                if ! is_integer "$mb"; then
                    echo "The value of -MB must be an integer."
                    usage
                fi
                shift
                shift
                ;;
            --debug)
                debug=true
                shift
                ;;
            --use-scr-need-checkpoint)
                use_scr_need_checkpoint=true
                shift
                ;;
            *)
                usage
                ;;
        esac
    done
fi

# Verify whether the required arguments were provided
if [ -z "$np" ] || [ -z "$nr" ] || [ -z "$nc" ] || [ -z "$qc" ]; then
    usage
fi

cd "$(dirname "$0")"
binary_filename=$(basename "$binary")

# Construct the execution command based on the arguments provided
command="mpirun -np $np"

if [ "$binary_filename" == "jacobi_noft" ]; then
    command+=" $binary_filename -p $nr -q $nc -NB $qc"
    if [ -n "$mb" ]; then
        command+=" -MB $mb"
    fi
elif [ "$binary_filename" == "jacobi_ulfm" ]; then
    command+=" --with-ft=ulfm --oversubscribe $binary_filename -p $nr -q $nc -NB $qc"
    if [ -n "$mb" ]; then
        command+=" -MB $mb"
    fi
else
    command+=" $binary_filename -p $nr -q $nc -NB $qc"
    if $debug; then
        command+=" --debug"
    fi
    if $use_scr_need_checkpoint; then
        command+=" --use-scr-need-checkpoint"
    fi
fi

# Execute the command constructed
eval "$command"
