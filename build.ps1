$ErrorActionPreference = 'Stop'

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location -LiteralPath $Root
try {
    & latexmk -xelatex -interaction=nonstopmode -halt-on-error -file-line-error main.tex
    if ($LASTEXITCODE -ne 0) {
        throw "latexmk failed with exit code $LASTEXITCODE"
    }

    $OutputDir = Join-Path $Root 'output\pdf'
    New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $Root 'main.pdf') `
        -Destination (Join-Path $OutputDir 'wigner_manybody_dynamics_draft.pdf') -Force
}
finally {
    Pop-Location
}

