#!/bin/bash

# Run script for CFD2D GUI application
# Quick launcher with common options

set -e

# Check if build exists
if [ ! -f "build/cfd2d_gui" ]; then
    echo "❌ Application not built. Run './build.sh' first."
    exit 1
fi

# Parse command line arguments
GUI_MODE=true
ARGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --no-gui)
            GUI_MODE=false
            ARGS="$ARGS $1"
            shift
            ;;
        --help)
            echo "CFD2D GUI - Quick Launcher"
            echo ""
            echo "Usage: ./run.sh [options]"
            echo ""
            echo "Options:"
            echo "  --no-gui              Run in headless mode"
            echo "  --Re=<value>          Set Reynolds number"
            echo "  --CFL=<value>         Set CFL condition"
            echo "  --solver=<pcg|mg>     Pressure solver type"
            echo "  --help                Show this help"
            echo ""
            echo "Examples:"
            echo "  ./run.sh                    # Run with GUI"
            echo "  ./run.sh --no-gui           # Run headless"
            echo "  ./run.sh --Re=100 --CFL=0.1 # Run with custom parameters"
            echo ""
            exit 0
            ;;
        *)
            ARGS="$ARGS $1"
            shift
            ;;
    esac
done

# Launch the application
if [ "$GUI_MODE" = true ]; then
    echo "🎮 Launching CFD2D GUI with dockable interface..."
    echo "   Use the top toolbar to toggle panels"
    echo "   Drag panels to rearrange the layout"
    echo "   Press 'Reset Layout' to restore default layout"
    echo ""
fi

./build/cfd2d_gui $ARGS 