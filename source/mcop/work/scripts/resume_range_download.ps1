param(
    [Parameter(Mandatory=$true)][string]$Url,
    [Parameter(Mandatory=$true)][string]$Path,
    [Parameter(Mandatory=$true)][Int64]$Start,
    [Parameter(Mandatory=$true)][Int64]$End
)

$parent = Split-Path -Parent $Path
New-Item -ItemType Directory -Force -Path $parent | Out-Null
$expected = $End - $Start + 1
$existing = if (Test-Path -LiteralPath $Path) { [Int64](Get-Item -LiteralPath $Path).Length } else { 0 }
if ($existing -gt $expected) { Remove-Item -LiteralPath $Path -Force; $existing = 0 }
$chunk = [Int64]2500000
while ($existing -lt $expected) {
    $requestStart = $Start + $existing
    $requestEnd = [Math]::Min($End, $requestStart + $chunk - 1)
    $tmp = "$Path.download"
    if (Test-Path -LiteralPath $tmp) { Remove-Item -LiteralPath $tmp -Force }
    curl.exe -L --fail --retry 4 --retry-all-errors --retry-delay 1 --max-time 55 -sS -A 'Mozilla/5.0' -r "$requestStart-$requestEnd" "$Url" -o "$tmp"
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    $got = [Int64](Get-Item -LiteralPath $tmp).Length
    $need = $requestEnd - $requestStart + 1
    if ($got -ne $need) { Write-Error "short range: got=$got need=$need start=$requestStart end=$requestEnd"; exit 2 }
    $outStream = [System.IO.File]::Open($Path, [System.IO.FileMode]::Append, [System.IO.FileAccess]::Write, [System.IO.FileShare]::Read)
    $inStream = [System.IO.File]::OpenRead($tmp)
    try { $inStream.CopyTo($outStream) } finally { $inStream.Dispose(); $outStream.Dispose() }
    Remove-Item -LiteralPath $tmp -Force
    $existing += $got
    Write-Output "$(Split-Path -Leaf $Path) $existing/$expected"
}
