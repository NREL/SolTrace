
# SolTrace Custom GPU Build Environment

This directory contains the configuration files required to generate the custom Docker container used for SolTrace's Continuous Integration (CI) pipeline on GPU-enabled runners.

## Purpose

Building SolTrace with GPU ray-tracing support (OptiX) requires the full NVIDIA CUDA Toolkit (several gigabytes) and the proprietary OptiX SDK. 

Instead of downloading and installing these heavy dependencies on the runner every time a test runs (which wastes CI minutes and budget), we pre-package them into a Docker image hosted on the GitHub Container Registry (GHCR):
`ghcr.io/nlr-soltrace/soltrace-build-env:latest`

When the GPU CI workflow (`CI.yml`) triggers, it pulls this image, providing a configured compilation environment to the runner hardware.

### Image
* **Base OS:** Ubuntu 24.04 LTS
* **CUDA Toolkit:** v13.1.0 (via official `nvidia/cuda` base image)
* **Build Tools:** `gcc`, `g++`, `cmake`, `git`
* **OptiX SDK:** v9.0.0 (Headers required for compilation)

---

## How to Update the Image

If you need to upgrade the CUDA version, install a newer OptiX SDK, or add new system dependencies, follow these steps to build and publish a new image.

### Step 1: Prepare the Files
1. Download the new Linux installer for the OptiX SDK (`.sh` file) from the NVIDIA Developer Portal.
2. Place the `.sh` file into this directory (`.github/docker/`).
3. Update the `Dockerfile` to reflect the new file names or base image versions.

### Step 2: Reactivate the Build Workflow
To prevent accidental builds and keep the repository clean, the build workflow is set to manual triggers only. You must temporarily re-enable the push trigger.
1. Open `.github/workflows/create-compiler-image.yml`.
2. Update the `on:` section to trigger on your current branch:
   ```yaml
   on:
     push:
       branches: [ your-active-branch ]
       paths:
         - '.github/docker/**'
   ```

### Step 3: Build and Publish
1. Commit the `.sh` file, your updated `Dockerfile`, and the modified `create-compiler-image.yml` file in a single commit.
2. Push to GitHub.
3. Navigate to the **Actions** tab on GitHub and watch the **Build and Publish compiler image** job execute. This job runs on a standard CPU runner and will publish the updated image to GHCR.

### Step 4: Clean Up
SDK installers are large binary files that bloat the repository history. Once the image is successfully published to GHCR, **delete the installer from the repository.**
1. Delete the raw OptiX `.sh` file from this directory.
2. Revert `.github/workflows/create-compiler-image.yml` to remove the `push` trigger (leaving only `workflow_dispatch`).
3. Commit and push these cleanup changes.
