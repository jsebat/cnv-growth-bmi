param([string]$src)
$pdf = $src -replace '\.docx$', '.pdf'
$w = New-Object -ComObject Word.Application; $w.Visible = $false
try {
  $doc = $w.Documents.Open($src, $false, $false)
  $n = 0; foreach ($p in $doc.Paragraphs) { $s = $p.Style.NameLocal; if ($s -like 'Heading*' -or $s -like 'heading*') { $p.Range.ListFormat.RemoveNumbers(); $n++ } }
  $t = 0; foreach ($tb in $doc.Tables) { $tb.Range.ParagraphFormat.Alignment = 0; $tb.AutoFitBehavior(0); $nc = $tb.Columns.Count; if ($nc -eq 2) { $tb.Columns.Item(1).Width = 130; $tb.Columns.Item(2).Width = 338 } elseif ($nc -eq 3) { $tb.Columns.Item(1).Width = 60; $tb.Columns.Item(2).Width = 230; $tb.Columns.Item(3).Width = 178 } else { $tb.AutoFitBehavior(2) }; $t++ }
  "headings un-numbered: $n; tables left-aligned and fitted: $t"
  $doc.Save(); $doc.ExportAsFixedFormat($pdf, 17); "pages: $($doc.ComputeStatistics(2))"; $doc.Close($false)
} finally { $w.Quit() }
$outDir = "C:\Users\JONATH~1\AppData\Local\Temp\claude\C--Users-Jonathan-Sebat-code-admixpgs\5d01a624-679d-4f8b-84b7-9e3205a6d729\scratchpad\pdfpages"
Remove-Item "$outDir\*.png" -Force -ErrorAction SilentlyContinue; $outDir = (Resolve-Path $outDir).Path
Add-Type -AssemblyName System.Runtime.WindowsRuntime
$null = [Windows.Data.Pdf.PdfDocument,Windows.Data.Pdf,ContentType=WindowsRuntime]; $null = [Windows.Storage.StorageFile,Windows.Storage,ContentType=WindowsRuntime]; $null = [Windows.Storage.StorageFolder,Windows.Storage,ContentType=WindowsRuntime]; $null = [Windows.Storage.Streams.IRandomAccessStream,Windows.Storage.Streams,ContentType=WindowsRuntime]
$asTaskGeneric = ([System.WindowsRuntimeSystemExtensions].GetMethods() | Where-Object { $_.Name -eq 'AsTask' -and $_.GetParameters().Count -eq 1 -and $_.GetParameters()[0].ParameterType.Name -eq 'IAsyncOperation`1' })[0]
$asTaskAction = ([System.WindowsRuntimeSystemExtensions].GetMethods() | Where-Object { $_.Name -eq 'AsTask' -and $_.GetParameters().Count -eq 1 -and $_.GetParameters()[0].ParameterType.Name -eq 'IAsyncAction' })[0]
function Await($op, $t) { $m = $asTaskGeneric.MakeGenericMethod($t); $task = $m.Invoke($null, @($op)); $task.Wait(-1) | Out-Null; $task.Result }
function AwaitAction($op) { $task = $asTaskAction.Invoke($null, @($op)); $task.Wait(-1) | Out-Null }
$file = Await ([Windows.Storage.StorageFile]::GetFileFromPathAsync($pdf)) ([Windows.Storage.StorageFile]); $pdfDoc = Await ([Windows.Data.Pdf.PdfDocument]::LoadFromFileAsync($file)) ([Windows.Data.Pdf.PdfDocument]); $folder = Await ([Windows.Storage.StorageFolder]::GetFolderFromPathAsync($outDir)) ([Windows.Storage.StorageFolder])
foreach ($i in 0..($pdfDoc.PageCount - 1)) { $page = $pdfDoc.GetPage([uint32]$i); $outFile = Await ($folder.CreateFileAsync("page-$($i+1).png", [Windows.Storage.CreationCollisionOption]::ReplaceExisting)) ([Windows.Storage.StorageFile]); $stream = Await ($outFile.OpenAsync([Windows.Storage.FileAccessMode]::ReadWrite)) ([Windows.Storage.Streams.IRandomAccessStream]); $opts = New-Object Windows.Data.Pdf.PdfPageRenderOptions; $opts.DestinationWidth = 900; AwaitAction ($page.RenderToStreamAsync($stream, $opts)); $stream.Dispose(); $page.Dispose() }
"rendered $($pdfDoc.PageCount) pages"
