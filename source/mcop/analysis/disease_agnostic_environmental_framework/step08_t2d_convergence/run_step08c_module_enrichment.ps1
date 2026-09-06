param(
    [string]$OutputDir = "$PSScriptRoot",
    [int]$SleepSeconds = 1
)

$ErrorActionPreference = "Stop"
$RawDir = Join-Path $OutputDir "raw_string"
$ModuleFile = Join-Path $OutputDir "t2d_step8c_network_modules.csv"
$StringBase = "https://string-db.org/api/json"

function Invoke-StringPost([hashtable]$Body, [string]$CachePath) {
    if (Test-Path $CachePath) {
        $existing = Get-Item $CachePath
        if ($existing.Length -gt 2) { return }
    }
    $lastError = $null
    for ($attempt = 1; $attempt -le 3; $attempt++) {
        try {
            $response = Invoke-WebRequest -Uri "$StringBase/enrichment" -Method Post -Body $Body -ContentType "application/x-www-form-urlencoded" -TimeoutSec 180 -UseBasicParsing
            [System.IO.File]::WriteAllText($CachePath, $response.Content, [System.Text.UTF8Encoding]::new($false))
            Start-Sleep -Seconds $SleepSeconds
            return
        } catch {
            $lastError = $_.Exception.Message
            if ($attempt -lt 3) { Start-Sleep -Seconds (5 * $attempt) }
        }
    }
    throw "STRING enrichment request failed after 3 attempts: $lastError"
}

$idToName = @{}
Get-ChildItem -LiteralPath $RawDir -File -Filter "*mapping.json" | ForEach-Object {
    $records = Get-Content $_.FullName -Raw | ConvertFrom-Json
    foreach ($record in $records) {
        if ($record.stringId -and $record.preferredName) {
            $idToName[[string]$record.stringId] = [string]$record.preferredName
        }
    }
}

$modules = Import-Csv $ModuleFile
$audit = @()
foreach ($module in $modules) {
    $ids = @($module.nodes -split ";" | Where-Object { $_ })
    $names = @($ids | ForEach-Object { if ($idToName.ContainsKey($_)) { $idToName[$_] } } | Sort-Object -Unique)
    $path = Join-Path $RawDir ("string_enrichment_{0}_{1}.json" -f $module.cluster_id, $module.module_id)
    if ($names.Count -ge 2) {
        Invoke-StringPost @{ identifiers = ($names -join "`n"); species = 9606 } $path
        $status = "queried"
    } else {
        "[]" | Set-Content -Encoding utf8 $path
        $status = "singleton_or_unmapped"
    }
    $audit += [ordered]@{ cluster_id = $module.cluster_id; module_id = $module.module_id; n_network_nodes = [int]$module.n_nodes; n_mapped_names = $names.Count; status = $status; cache = $path }
}
$audit | ConvertTo-Json -Depth 8 | Set-Content -Encoding utf8 (Join-Path $OutputDir "STEP8C_MODULE_ENRICHMENT_AUDIT.json")
Write-Output ("Module enrichment complete: {0} modules" -f $audit.Count)
