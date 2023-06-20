#!/bin/sh
network_3D_subnet_key_value="3D_subnet"

network_3D_distance_key_value="3D_distance"
if [ $# -lt 1 -o $# -gt 10 ]; then
        echo 1>&2 "Usage: run.sh --network_config="" --mri_network(OP)="" --3D_subnet(OP)="" --3D_distance(OP)="" "
    exit 127
fi

subdir="$PWD"

network_config_file_key_value="network_config"
params_file_key_value="params_config"
output_location="output"
network_mri_file_key_value="mri_network"
network_3D_subnet_key_value="3D_subnet"

network_3D_distance_key_value="3D_distance"

mri_run=""
3D_run=""

while test $# -gt 0;do
    case "$1" in
       -$network_config_file_key_value)
               shift
               network_config_file=$1
               network_run=true
               shift
               ;;
       -$params_file_key_value)
               shift
               params_config_value=$1
               shift
               ;;
       -$output_location)
               shift
               output_location_value=$1
               shift
               ;;
       -$network_mri_file_key_value)
               shift
               network_mri_file_value=$1
               mri_run="true"
               shift
               ;;
                      -$network_3D_subnet_key_value)
               shift
               network_3D_subnet_value=$1
               3D_run="true"
               shift
               ;;
       -$network_3D_distance_key_value)
               shift
               network_3D_distance_value=$1
               shift
               ;;
       *)
               return 1;
               ;;
    esac
done

#wrk_dir="/home/scigap/NeuroTest/workdir"
binary_location="/home/scigap/NeuroTest/binaries/"

network_binary="$binary_location""generate_network"
mindcurrent_binary="$binary_location""mindcurrent"

if [[ $mri_run = true ]]; then
  echo "MRI and Network."
 # cd $wrk_dir
  $network_binary  $network_config_file $network_mri_file_value > connection_info2
elif [[ 3D_run = true ]]; then
  echo "Network and 3D"
  #cd $wrk_dir
  $network_binary $network_config_file $network_3D_subnet_value $network_3D_distance_key_value > connection_info2
elif [[ $network_run = true ]]; then
  echo "Network only"
  #cd $wrk_dir
  $network_binary $network_config_file > connection_info2
else
  echo "Invalid Input paramerters"
fi

$mindcurrent_binary $params_config_value $output_location_value connection_info2
