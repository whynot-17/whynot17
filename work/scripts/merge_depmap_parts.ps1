$dir = Join-Path (Get-Location) 'work\phase7b_depmap\raw'
function Merge-Parts([string]$prefix, [string]$out) {
    $target = Join-Path $dir $out
    if (Test-Path -LiteralPath $target) { Remove-Item -LiteralPath $target -Force }
    $fs = [IO.File]::OpenWrite($target)
    try {
        Get-ChildItem -LiteralPath $dir -Filter "$prefix.part??" | Sort-Object Name | ForEach-Object {
            $inStream = [IO.File]::OpenRead($_.FullName)
            try { $inStream.CopyTo($fs) } finally { $inStream.Dispose() }
        }
    } finally { $fs.Dispose() }
    Write-Output "$out $((Get-Item -LiteralPath $target).Length)"
}
Merge-Parts 'crispr' 'CRISPRGeneEffect.csv'
Merge-Parts 'expr' 'OmicsExpressionProteinCodingGenesTPMLogp1.csv'
