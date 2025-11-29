# MATLAB + Communications Toolbox Support Package for USRP Installation Guide

**Platforms:** Windows, Linux, macOS

**Time Estimate:** 30–60 minutes

---

## Prerequisites

- **MATLAB Version:** R2020a or later
- **Communications Toolbox:** Must be installed
- **RAM:** 8 GB minimum (16 GB recommended)
- **Disk Space:** 10 GB free (for MATLAB + support package)
- **USRP Device:** B200, B210, N210, X310, or similar
- **Internet:** Required for package downloads

---

## Step 1: Verify MATLAB and Communications Toolbox

Open MATLAB and verify that Communications Toolbox is installed:

```matlab
ver  % Check installed products
```

Look for "Communications Toolbox" in the list. If not installed, purchase and install it via MATLAB's Add-Ons manager.

---

## Step 2: Install the Support Package for USRP Radio

### Method A: Using the Add-Ons Manager (Recommended)

1. Open MATLAB.
2. Click **Home** → **Add-Ons** → **Get Add-Ons**.
3. Search for `"Communications Toolbox Support Package for USRP Radio"`.
4. Click **Install** and follow the prompts.

The installer will automatically download and configure the support package.

### Method B: Using Command Line

In the MATLAB command window:

```matlab
supportPackageInstaller
```

Select "Communications Toolbox Support Package for USRP Radio" from the list and install.

---

## Step 3: Install UHD (USRP Hardware Driver)

UHD is required for MATLAB to communicate with USRP devices.

### Windows

1. Download the latest UHD binary from:  
   **https://files.ettus.com/binaries/uhd/**

2. Choose the Windows 64-bit installer (e.g., `uhd_4.0.0.0_Win64.exe`).

3. Run the installer as **Administrator** and follow the setup wizard.

4. Add UHD to your system PATH (usually done automatically):
   - Typical path: `C:\Program Files\UHD\bin`
   - If not added automatically:
     - Right-click **This PC** → **Properties** → **Advanced System Settings** → **Environment Variables**.
     - Add the UHD `bin` directory to the `PATH` variable.
     - Restart MATLAB.

### Linux (Ubuntu / Debian)

If you did not already install UHD from the GNU Radio guide, install it via apt:

```bash
sudo apt-get install libuhd-dev uhd-host
uhd_images_downloader
```

Set the library path (add to `~/.bashrc`):

```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

Then reload:

```bash
source ~/.bashrc
```

### macOS

Install UHD via Homebrew:

```bash
brew install uhd
uhd_images_downloader
```

---

## Step 4: Download FPGA Images

FPGA images are required for your USRP device to function. Run this command in a terminal:

```bash
uhd_images_downloader
```

This downloads the necessary firmware for your specific USRP model.

---

## Step 5: Verify UHD Installation

From a terminal (Windows PowerShell, Linux bash, or macOS terminal), verify UHD:

**Windows (PowerShell):**

```powershell
uhd_find_devices
uhd_config_info --version
```

**Linux / macOS:**

```bash
uhd_find_devices
uhd_config_info --version
```

**Expected output:**
```
UHD 4.0.0.0 or similar version
```

---

## Step 6: Test MATLAB-USRP Connection

Connect your USRP device via USB 3.0 (or Ethernet for supported models) and verify MATLAB can detect it.

Create a test script (e.g., `test_connection.m`) and run it in MATLAB:

```matlab
%% Test USRP Connection
clear; clc; close all;

try
    % Create USRP receiver object
    radio = comm.SDRuReceiver('Platform', 'B210');
    
    % Display information
    fprintf('USRP Device Connected Successfully!\n');
    fprintf('Platform: %s\n', radio.Platform);
    
    % Get device info
    devInfo = info(radio);
    fprintf('Status: %s\n', devInfo.Status);
    
    % Release resources
    release(radio);
    
catch ME
    fprintf('Error: %s\n', ME.message);
    fprintf('Troubleshooting:\n');
    fprintf('  - Check USB 3.0 connection\n');
    fprintf('  - Verify UHD drivers are installed\n');
    fprintf('  - Run uhd_find_devices from terminal\n');
end
```

**Expected output:**
```
USRP Device Connected Successfully!
Platform: B210
Status: OK
```

---

## Troubleshooting

### "Device Not Found" Error

1. Check USB 3.0 cable and port (USB 2.0 may not work reliably).
2. Power on the USRP device.
3. Verify UHD installation:
   ```bash
   uhd_find_devices
   ```
4. Download FPGA images:
   ```bash
   uhd_images_downloader
   ```

### "Permission Denied" (Linux)

Add your user to the `usrp` group:

```bash
sudo usermod -a -G usrp $USER
newgrp usrp
```

Restart MATLAB and reconnect the device.

### UHD Not Found in PATH (Windows)

1. Verify UHD is installed in `C:\Program Files\UHD\` (or custom location).
2. Add `C:\Program Files\UHD\bin` to Windows PATH:
   - **Settings** → **System** → **About** → **Advanced System Settings**.
   - **Environment Variables** → **PATH** → **Edit** → **New**.
   - Add `C:\Program Files\UHD\bin`.
   - Restart MATLAB.

### MATLAB Support Package Installation Fails

- Ensure you have at least 2 GB free disk space.
- Check your internet connection.
- Try restarting MATLAB and attempting installation again.
- If using a proxy, configure MATLAB's proxy settings in Preferences.

---

## Next Steps

- Explore the MATLAB **Communications Toolbox** documentation for USRP examples.
- Review sample scripts in your workspace or MATLAB documentation for:
  - Signal reception
  - Signal transmission
  - Loopback testing
  - Spectrum analysis
- Test with simple signal reception before attempting transmission.

---

## Verification Checklist

- [ ] MATLAB and Communications Toolbox installed
- [ ] Support Package for USRP installed (Add-Ons)
- [ ] UHD drivers installed (Windows / Linux / macOS)
- [ ] FPGA images downloaded (`uhd_images_downloader`)
- [ ] USRP device connected and powered
- [ ] `uhd_find_devices` detects the device
- [ ] MATLAB test script runs without errors

---

## Additional Resources

- **MATLAB USRP Documentation:** https://www.mathworks.com/help/comm/
- **UHD Manual:** https://files.ettus.com/manual/
- **Ettus Research USRP:** https://www.ettus.com/
- **MATLAB Answers:** https://www.mathworks.com/matlabcentral/answers/

---

**Last Updated:** November 2025
