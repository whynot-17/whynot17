param(
    [string]$Message = "Update analysis outputs"
)

$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Set-Location -LiteralPath $repoRoot

$branch = (git branch --show-current).Trim()
if ($branch -ne "main") {
    throw "Expected branch 'main', found '$branch'."
}

git add .
git diff --cached --quiet
if ($LASTEXITCODE -eq 0) {
    Write-Host "No new tracked changes to push."
    exit 0
}

git commit -m $Message
git push origin main
Write-Host "Pushed to origin/main."
