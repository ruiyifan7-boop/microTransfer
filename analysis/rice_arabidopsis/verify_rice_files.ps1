$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$DataDir = Join-Path $Root "data\rice"
$Manifest = Import-Csv -Delimiter "`t" (Join-Path $DataDir "edwards_rice_dryad_manifest.tsv")

$rows = foreach ($row in $Manifest) {
    $path = Join-Path $DataDir $row.file
    if (-not (Test-Path -LiteralPath $path)) {
        [PSCustomObject]@{
            File = $row.file
            Status = "MISSING"
            Bytes = ""
            MD5 = ""
        }
        continue
    }

    $item = Get-Item -LiteralPath $path
    $hash = (Get-FileHash -Algorithm MD5 -LiteralPath $path).Hash.ToLowerInvariant()
    $ok = ($item.Length -eq [int64]$row.expected_bytes) -and ($hash -eq $row.md5)
    [PSCustomObject]@{
        File = $row.file
        Status = if ($ok) { "OK" } else { "FAILED" }
        Bytes = $item.Length
        MD5 = $hash
    }
}

$rows | Format-Table -AutoSize
if ($rows.Status -contains "MISSING" -or $rows.Status -contains "FAILED") {
    exit 2
}

