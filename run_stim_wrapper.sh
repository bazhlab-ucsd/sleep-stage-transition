#!/bin/sh

# Script Configuration
binary_location="/path/to/binaries"
network_binary="$binary_location/generate_network"
mindcurrent_binary="$binary_location/mindcurrent"
output_location="output"

# Key Values
network_config_file_key="network_config"
params_file_key="params_config"
mri_network_key="mri_network"
subnet_key="3D_subnet"
distance_key="3D_distance"

# Variables
network_config_file=""
params_file=""
output_dir=""
mri_file=""
subnet_file=""
distance_file=""
network_run=false
mri_run=false
subnet_run=false

# Usage Information
usage() {
  echo "Usage: $0"
  echo "  --network_config <file>    Required: Network configuration file"
  echo "  --params_config <file>     Required: Parameters file"
  echo "  --output <dir>             Optional: Output directory (default: output)"
  echo "  --mri_network <file>       Optional: MRI network file"
  echo "  --3D_subnet <file>         Optional: 3D subnet file"
  echo "  --3D_distance <file>       Optional: 3D distance file"
  exit 1
}

# Parse Command-Line Arguments
while [ $# -gt 0 ]; do
  case "$1" in
    --network_config)
      network_config_file="$2"
      network_run=true
      shift 2
      ;;
    --params_config)
      params_file="$2"
      shift 2
      ;;
    --output)
      output_dir="$2"
      shift 2
      ;;
    --mri_network)
      mri_file="$2"
      mri_run=true
      shift 2
      ;;
    --3D_subnet)
      subnet_file="$2"
      subnet_run=true
      shift 2
      ;;
    --3D_distance)
      distance_file="$2"
      shift 2
      ;;
    *)
      echo "Invalid argument: $1"
      usage
      ;;
  esac
done

# Check Required Arguments
if [ -z "$network_config_file" ] || [ -z "$params_file" ]; then
  echo "Error: --network_config and --params_config are required."
  usage
fi

# Set Default Output Directory
output_dir="${output_dir:-$output_location}"

# Verify File Existence
if [ ! -f "$network_config_file" ]; then
  echo "Error: Network config file not found: $network_config_file"
  exit 1
fi

if [ ! -f "$params_file" ]; then
  echo "Error: Params config file not found: $params_file"
  exit 1
fi

# Execute the Appropriate Workflow
connection_info="$output_dir/connection_info2"
mkdir -p "$output_dir"

if [ "$mri_run" = true ]; then
  echo "Running MRI and Network workflow..."
  "$network_binary" "$network_config_file" "$mri_file" > "$connection_info"
elif [ "$subnet_run" = true ]; then
  echo "Running Network and 3D workflow..."
  "$network_binary" "$network_config_file" "$subnet_file" "$distance_file" > "$connection_info"
elif [ "$network_run" = true ]; then
  echo "Running Network-only workflow..."
  "$network_binary" "$network_config_file" > "$connection_info"
else
  echo "Error: Invalid combination of parameters."
  usage
fi

# Run the Mindcurrent Binary
"$mindcurrent_binary" "$params_file" "$output_dir" "$connection_info"

echo "Execution completed successfully."