#!/bin/bash

################################################################################
# Script: execute.sh
# Description: This script executes the binary with the provided arguments.
#
# Usage: execute.sh <binary> < -f <flag_file> | -np <NP> -p <NR> -q <NC>
#           -NB <QC> [-MB <QR>] [--hostfile=<hostfile>] [--debug]
#           [--use-scr-need-checkpoint] [--run-until-success]
#           [--retry-delay <delay_time>] >
#
#   <binary>    Binary to be executed
#   <flag_file> File with the flags to be used by the application
#   <NP>        Number of processes used by the application
#   <NR>        Number of partitions along the x-axis.
#   <NC>        Number of partitions along the y-axis.
#   <QC>        Number of columns in the input matrix.
#   <QR>        Number of rows in the input matrix.
#   <hostfile>  Hostfile for mpirun
#   --debug                    (only for jacobi_scr) Print debug information
#   --use-scr-need-checkpoint  (only for jacobi_scr) Use SCR_Need_checkpoint
#   --run-until-success        (only for jacobi_scr) If specified, run the
#       command until it returns 0  (caution: this may cause an infinite loop)
#   --retry-delay <delay_time> (only if --run-until-success is specified) Delay
#       time in seconds between retries. Default: 1 second.
#
# Usage examples:
#   ./execute.sh jacobi_noft -f flags.conf
#   ./execute.sh jacobi_ulfm -np 4 -p 2 -q 2 -NB 1000 -MB 1000 --hostfile hostfile
#   ./execute.sh jacobi_scr -np 4 -p 2 -q 2 -NB 1000 --hostfile hostfile
#       --debug --use-scr-need-checkpoint --run-until-success --retry-delay 0.5
#
# Notes:
#   - When jacobi_ulfm is used, the --oversubscribe and --with-ft=ulfm flags are
#   added to the mpirun command.
#
################################################################################

# print usage
usage() {
    echo "Usage: $0 <binary> <(-f <flag_file>) | (-np <NP> -p <NR> -q <NC> -NB <QC> [-MB <QR>] [--hostfile=<hostfile>] [--debug] [--use-scr-need-checkpoint] [--run-until-success] [--retry-delay <delay_time>])>"
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
            qr=$value
            if ! is_integer "$qr"; then
                echo "The value of -MB must be an integer."
                usage
            fi
            ;;
        --hostfile=*)
            hostfile="${flag#*=}"
            ;;
        --debug)
            debug=true
            ;;
        --use-scr-need-checkpoint)
            use_scr_need_checkpoint=true
            ;;
        --run-until-success)
            run_until_success=true
            ;;
        --retry-delay)
            delay_time=$value
            ;;
        *)
            unknown_flags+=("$flag")
            ;;
        esac
    done <"$1"

    # Verify whether there are unknown flags and print them
    if [ ${#unknown_flags[@]} -gt 0 ]; then
        for flag in "${unknown_flags[@]}"; do
            echo "Unknown flag: $flag"
        done
        usage
    fi
}

# Construct the command to be executed
construct_command() {
    # mpirun config
    command="mpirun -np $np"

    if [ -n "$hostfile" ]; then
        command+=" --hostfile=$hostfile"
    fi

    if [ -n "$extra_args" ]; then
        command+=" $extra_args"
    fi

    # binary config
    command+=" $binary -p $nr -q $nc -NB $qc"

    if [ -n "$qr" ]; then
        command+=" -MB $qr"
    fi

    if [ "$debug" = "true" ]; then
        command+=" --debug"
    fi

    if [ "$use_scr_need_checkpoint" = "true" ]; then
        command+=" --use-scr-need-checkpoint"
    fi
}

# Execute the command until it returns 0
run_command() {
    if [ "$run_until_success" = "false" ]; then
        $command
    else
        while true; do
            $command
            if [ $? -eq 0 ]; then
                break
            fi
            sleep "$delay_time"
        done
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
run_until_success=false
delay_time=1
command=""

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
            qr=$2
            if ! is_integer "$qr"; then
                echo "The value of -MB must be an integer."
                usage
            fi
            shift
            shift
            ;;
        --hostfile=*)
            hostfile="${1#*=}"
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
        --run-until-success)
            run_until_success=true
            shift
            ;;
        --retry-delay)
            echo "delay_time: $2"
            delay_time=$2
            shift
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

if [ "$binary_filename" == "jacobi_ulfm" ]; then
    extra_args="--with-ft=ulfm --oversubscribe"
else
    extra_args=""
fi

construct_command
run_command
