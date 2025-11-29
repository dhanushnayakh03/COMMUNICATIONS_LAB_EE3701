# GNU Radio 3.8 + UHD 4.0 + gr-ettus Installation Guide

**Platform:** Ubuntu 20.04 LTS

**Time Estimate:** 2–3 hours (including build time)

---

## Prerequisites

- **OS:** Ubuntu 20.04 LTS
- **RAM:** 8 GB minimum (16 GB recommended)
- **Disk Space:** 50 GB free
- **Internet:** Required for package downloads and cloning repositories
- **USRP Device:** B200, B210, N210, X310, or similar

---

## Step 1: System Update

Update your package manager and install runtime libraries for building:

```bash
sudo apt-get update && sudo apt-get upgrade -y
```

---

## Step 2: Install Build Dependencies

Install all required libraries and compilers for UHD, GNU Radio, and gr-ettus:

```bash
sudo apt install -y git cmake g++ libboost-all-dev libgmp-dev \
    swig python3-numpy python3-mako python3-sphinx python3-lxml \
    doxygen libfftw3-dev libsdl1.2-dev libgsl-dev libqwt-qt5-dev \
    libqt5opengl5-dev python3-pyqt5 liblog4cpp5-dev libzmq3-dev \
    python3-yaml python3-click python3-click-plugins python3-zmq \
    python3-scipy python3-gi python3-gi-cairo gobject-introspection \
    gir1.2-gtk-3.0 build-essential libusb-1.0-0-dev python3-docutils \
    python3-setuptools python3-ruamel.yaml python-is-python3 \
    libtinfo5 libncurses5
```

---

## Step 3: (Optional) Xilinx Vivado for RFNOC

If you plan to use RFNOC (FPGA development), prepare the Vivado installation directory:

```bash
sudo mkdir -p /opt/Xilinx
sudo chmod -R 777 /opt/Xilinx
```

Download Xilinx Vivado 2019.1 from the Xilinx website and run the installer:

```bash
cd ~/Downloads/Xilinx_Vivado_SDK_2019.1_0524_1430/
./xsetup
```

Point the installer to `/opt/Xilinx` and select **Vivado HL Design Edition**. Installation may take 1–2 hours.

---

## Step 4: Build and Install UHD 4.0

Clone the UHD repository and build it from source:

```bash
git clone --branch UHD-4.0 https://github.com/ettusresearch/uhd.git uhd
cd uhd/host
mkdir -p build && cd build
cmake ..
make -j$(nproc)
make test
sudo make install
sudo ldconfig
```

**Note:** The `make` step may take 30–60 minutes depending on your system. `-j$(nproc)` uses all available CPU cores.

---

## Step 5: Build and Install GNU Radio 3.8

Clone and build GNU Radio from the maintenance branch:

```bash
git clone --branch maint-3.8 --recursive https://github.com/gnuradio/gnuradio.git gnuradio
cd gnuradio
mkdir -p build && cd build
cmake ..
make -j$(nproc)
make test
sudo make install
sudo ldconfig
```

---

## Step 6: Build and Install gr-ettus

Build gr-ettus with Qt support for gr-fosphor (spectrum visualization):

```bash
git clone --branch maint-3.8-uhd4.0 https://github.com/ettusresearch/gr-ettus.git gr-ettus
cd gr-ettus
mkdir -p build && cd build
cmake -DENABLE_QT=True ..
make -j$(nproc)
make test
sudo make install
```

---

## Step 7: Configure Python Path

Ensure Python can find the installed modules. Add the following to your `~/.bashrc`:

```bash
echo 'export PYTHONPATH=/usr/local/lib/python3/dist-packages/:$PYTHONPATH' >> ~/.bashrc
source ~/.bashrc
```

Also run:

```bash
sudo ldconfig
```

---

## Step 8: Verify Installation

Check that UHD, GNU Radio, and gr-ettus are installed correctly:

```bash
# Check UHD version
uhd_config_info --version

# Expected output: UHD 4.0.0.0 or similar

# Check GNU Radio version
gnuradio-config-info --version

# Expected output: v3.8.x.x or similar

# Check for RFNOC support
rfnocmodtool help
```

---

## Step 9: Download FPGA Images

Download FPGA images for your USRP device. This is required for the device to function:

```bash
uhd_images_downloader
```

---

## Step 10: Test Device Connection

Connect your USRP via USB 3.0 or Ethernet and run:

```bash
uhd_find_devices
```

**Expected output:** You should see your USRP device listed with its serial number and type, e.g.:
```
--------------------------------------------------
-- UHD Device 0
--------------------------------------------------
Device Address:
    type: b210
    serial: XXXXXX
```

---

## Troubleshooting

### Permission Denied on Device

Add your user to the `usrp` group and restart:

```bash
sudo usermod -a -G usrp $USER
newgrp usrp
```

Then reconnect the USRP device or restart MATLAB/MATLAB GUI.

### Build Fails

- Check `CMakeError.log` in the `build/` directory for missing libraries.
- Ensure all `-dev` packages from Step 2 are installed.
- Try removing the build directory and reconfiguring:

```bash
rm -rf build
mkdir build && cd build
cmake ..
```

### Device Not Found

- Verify USB 3.0 connection (USRP may not work reliably on USB 2.0).
- Check that the device is powered on.
- Run `uhd_images_downloader` to update FPGA images.

### Python Module Import Errors

Ensure `PYTHONPATH` is set (Step 7) and run:

```bash
python3 -c "import uhd; print(uhd.__version__)"
```

---

## Next Steps

- Explore GNU Radio Companion (GRC) GUI: `gnuradio-companion`
- Check out sample applications in the GNU Radio documentation.
- Run `usrp_test.py`


---

## Additional Resources

- **UHD Documentation:** https://files.ettus.com/manual/
- **GNU Radio Documentation:** https://www.gnuradio.org/
- **gr-ettus Repository:** https://github.com/ettusresearch/gr-ettus

---

**Last Updated:** November 2025
