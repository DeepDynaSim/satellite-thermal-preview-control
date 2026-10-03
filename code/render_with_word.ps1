param([Parameter(Mandatory=$true)][string]$DocxPath,[Parameter(Mandatory=$true)][string]$PdfPath,[string]$AuditPath)
$ErrorActionPreference='Stop'
$resolvedDocx=(Resolve-Path -LiteralPath $DocxPath).Path
$resolvedPdf=[System.IO.Path]::GetFullPath($PdfPath)
$wordApp=$null
$wordDoc=$null
try {
    $wordApp=New-Object -ComObject Word.Application
    $wordApp.Visible=$false
    $wordApp.DisplayAlerts=0
    $wordDoc=$wordApp.Documents.Open($resolvedDocx,$false,$false)
    for ($i=1;$i -le $wordDoc.OMaths.Count;$i++) {
        $equation=$wordDoc.OMaths.Item($i)
        $equation.Range.Font.Name='Cambria Math'
        $equation.BuildUp()
    }
    $wordDoc.Fields.Update() | Out-Null
    $wordDoc.Repaginate()
    $wordDoc.Save()
    $audit=[ordered]@{word_version=$wordApp.Version;page_count=$wordDoc.ComputeStatistics(2);word_count=$wordDoc.ComputeStatistics(0);native_OMath_count=$wordDoc.OMaths.Count;table_count=$wordDoc.Tables.Count;inline_shapes=$wordDoc.InlineShapes.Count}
    $wordDoc.ExportAsFixedFormat($resolvedPdf,17)
    if ($AuditPath) { $audit | ConvertTo-Json | Set-Content -LiteralPath $AuditPath -Encoding utf8 }
    $audit | ConvertTo-Json
} finally {
    if ($wordDoc) { $wordDoc.Close(0); [System.Runtime.InteropServices.Marshal]::ReleaseComObject($wordDoc) | Out-Null }
    if ($wordApp) { $wordApp.Quit(); [System.Runtime.InteropServices.Marshal]::ReleaseComObject($wordApp) | Out-Null }
}
