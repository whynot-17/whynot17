param(
    [string]$OutputDir = "$PSScriptRoot",
    [int]$RequiredScore = 700,
    [string]$NetworkType = "functional",
    [int]$MaxPostGenes = 2000,
    [int]$SleepSeconds = 1
)

$ErrorActionPreference = "Stop"
$Step7Membership = Join-Path (Split-Path $OutputDir -Parent) "step07_genecard_convergence\t2d_cluster_ctd_gene_membership.csv"
$Step7Membership = [System.IO.Path]::GetFullPath($Step7Membership)
$RawDir = Join-Path $OutputDir "raw_string"
New-Item -ItemType Directory -Force -Path $RawDir | Out-Null

$TierA = @("cluster_11", "cluster_5", "cluster_6", "cluster_8")
$StringBase = "https://string-db.org/api/json"

function Get-Genes([string]$ClusterId) {
    $rows = Import-Csv $Step7Membership
    return @($rows | Where-Object { $_.cluster_id -eq $ClusterId } | ForEach-Object { $_.gene_symbol } | Sort-Object -Unique)
}

function Invoke-StringPost([string]$Endpoint, [hashtable]$Body, [string]$CachePath) {
    if (Test-Path $CachePath) {
        $existing = Get-Item $CachePath
        if ($existing.Length -gt 2) { return }
    }
    $lastError = $null
    for ($attempt = 1; $attempt -le 3; $attempt++) {
        try {
            $response = Invoke-WebRequest -Uri "$StringBase/$Endpoint" -Method Post -Body $Body -ContentType "application/x-www-form-urlencoded" -TimeoutSec 180 -UseBasicParsing
            [System.IO.File]::WriteAllText($CachePath, $response.Content, [System.Text.UTF8Encoding]::new($false))
            Start-Sleep -Seconds $SleepSeconds
            return
        } catch {
            $lastError = $_.Exception.Message
            if ($attempt -lt 3) { Start-Sleep -Seconds (5 * $attempt) }
        }
    }
    throw "STRING request failed after 3 attempts: $Endpoint; $lastError"
}

function Get-Blocks([string[]]$Genes) {
    if ($Genes.Count -le $MaxPostGenes) {
        return [pscustomobject]@{ block_index = 0; genes = [string[]]$Genes }
    }
    $nBlocks = [math]::Ceiling($Genes.Count / 1000)
    $blockSize = [math]::Ceiling($Genes.Count / $nBlocks)
    $blocks = @()
    for ($i = 0; $i -lt $Genes.Count; $i += $blockSize) {
        $end = [math]::Min($i + $blockSize, $Genes.Count)
        $blocks += [pscustomobject]@{ block_index = $blocks.Count; genes = [string[]]($Genes[$i..($end - 1)]) }
    }
    return $blocks
}

function Get-Union([string[]]$Left, [string[]]$Right) {
    return @($Left + $Right | Sort-Object -Unique)
}

$manifest = [ordered]@{
    status = "acquisition_started"
    endpoint_base = $StringBase
    species = 9606
    required_score = $RequiredScore
    network_type = $NetworkType
    add_nodes = 0
    input_source = $Step7Membership
    input_rule = "Step 7 overlap genes only for Tier A axes; frozen 11-cluster union used only as empirical randomization background"
    generated_at = (Get-Date).ToUniversalTime().ToString("o")
    axes = @{}
    background = @{}
}

foreach ($cluster in $TierA) {
    $genes = Get-Genes $cluster
    $mappingPath = Join-Path $RawDir "string_${cluster}_mapping.json"
    Invoke-StringPost "get_string_ids" @{ identifiers = ($genes -join "`n"); species = 9606 } $mappingPath
    $blocks = @(Get-Blocks $genes)
    $requests = @()
    if ($blocks.Count -eq 1) {
        $path = Join-Path $RawDir "string_${cluster}_network_b001.json"
        Invoke-StringPost "network" @{ identifiers = ($blocks[0].genes -join "`n"); species = 9606; required_score = $RequiredScore; network_type = $NetworkType; add_nodes = 0 } $path
        $requests += [ordered]@{ id = "b001"; n_genes = $blocks[0].genes.Count; block_a = 0; block_b = 0; path = $path }
    } else {
        $requestIndex = 0
        for ($i = 0; $i -lt $blocks.Count; $i++) {
            for ($j = $i + 1; $j -lt $blocks.Count; $j++) {
                $requestIndex++
                $union = Get-Union $blocks[$i].genes $blocks[$j].genes
                $requestId = "b{0:D3}" -f $requestIndex
                $path = Join-Path $RawDir "string_${cluster}_network_${requestId}.json"
                Invoke-StringPost "network" @{ identifiers = ($union -join "`n"); species = 9606; required_score = $RequiredScore; network_type = $NetworkType; add_nodes = 0 } $path
                $requests += [ordered]@{ id = $requestId; n_genes = $union.Count; block_a = $i; block_b = $j; path = $path }
            }
        }
    }
    $manifest.axes[$cluster] = [ordered]@{ n_input_genes = $genes.Count; n_blocks = $blocks.Count; n_requests = $requests.Count; mapping = $mappingPath; requests = $requests }
}

$allRows = Import-Csv $Step7Membership
$backgroundGenes = @($allRows | ForEach-Object { $_.gene_symbol } | Sort-Object -Unique)
$backgroundMappingPath = Join-Path $RawDir "string_background_11clusters_mapping.json"
Invoke-StringPost "get_string_ids" @{ identifiers = ($backgroundGenes -join "`n"); species = 9606 } $backgroundMappingPath
$backgroundBlocks = @(Get-Blocks $backgroundGenes)
$backgroundRequests = @()
$requestIndex = 0
if ($backgroundBlocks.Count -eq 1) {
    $path = Join-Path $RawDir "string_background_11clusters_network_b001.json"
    Invoke-StringPost "network" @{ identifiers = ($backgroundBlocks[0].genes -join "`n"); species = 9606; required_score = $RequiredScore; network_type = $NetworkType; add_nodes = 0 } $path
    $backgroundRequests += [ordered]@{ id = "b001"; n_genes = $backgroundBlocks[0].genes.Count; block_a = 0; block_b = 0; path = $path }
} else {
    for ($i = 0; $i -lt $backgroundBlocks.Count; $i++) {
        for ($j = $i + 1; $j -lt $backgroundBlocks.Count; $j++) {
            $requestIndex++
            $union = Get-Union $backgroundBlocks[$i].genes $backgroundBlocks[$j].genes
            $requestId = "b{0:D3}" -f $requestIndex
            $path = Join-Path $RawDir "string_background_11clusters_network_${requestId}.json"
            Invoke-StringPost "network" @{ identifiers = ($union -join "`n"); species = 9606; required_score = $RequiredScore; network_type = $NetworkType; add_nodes = 0 } $path
            $backgroundRequests += [ordered]@{ id = $requestId; n_genes = $union.Count; block_a = $i; block_b = $j; path = $path }
        }
    }
}
$manifest.background = [ordered]@{ n_input_genes = $backgroundGenes.Count; n_blocks = $backgroundBlocks.Count; n_requests = $backgroundRequests.Count; mapping = $backgroundMappingPath; requests = $backgroundRequests }
$manifest.status = "acquisition_complete"
$manifestPath = Join-Path $OutputDir "STEP8C_STRING_ACQUISITION_MANIFEST.json"
$manifest | ConvertTo-Json -Depth 12 | Set-Content -Encoding utf8 $manifestPath
Write-Output ("STRING acquisition complete: {0}" -f $manifestPath)
