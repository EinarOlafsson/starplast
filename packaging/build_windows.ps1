# Build starplast.exe.  powershell -File packaging\build_windows.ps1
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
$Version = python -c "import tomllib,pathlib;print(tomllib.loads(pathlib.Path(r'$Root/pyproject.toml').read_text())['project']['version'])"
Remove-Item -Recurse -Force "$Root\build\win" -ErrorAction SilentlyContinue
pyinstaller --noconfirm --distpath "$Root\build\win\dist" --workpath "$Root\build\win\work" `
            "$Root\packaging\starplast.spec"
Write-Host "built: build\win\dist\starplast\starplast.exe (version $Version)"
# Inno Setup, when present, wraps it into a single-file installer.
$iscc = "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe"
if (Test-Path $iscc) { & $iscc "$Root\packaging\starplast.iss" }
