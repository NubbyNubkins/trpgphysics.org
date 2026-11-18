# trpgphysics.org

The static website for Dynamic Relational Transitioning Geometry (TRPG): https://trpgphysics.org

Everything under `site/` is the published site, file for file. The research code, data and paper archive live in
the companion repository [TRPG-Unified-Framework](https://github.com/NubbyNubkins/TRPG-Unified-Framework).

## How a change goes live

A change merged to `main` that touches `site/` is picked up by the web host within about ten minutes. The host
pulls from this repository; nothing is pushed to it. Before a new version is published the host:

* looks only inside `site/`, and refuses symlinks, hidden files, unexpected file types and oversized trees;
* checks every line of `site/downloads/SHA256SUMS.txt` against the files;
* switches over only if the key pages answer correctly, and switches back by itself if they do not.

## Editing

Edit files under `site/`. If you change a file listed in `site/downloads/SHA256SUMS.txt`, update its line there too,
or the update is refused. Static files only (html, css, js, json, txt, xml, csv, tsv, py, zip, pdf, docx, xlsx,
png, svg and similar).

## Licence

See `LICENSE` (the same text as `site/LICENSE.txt`).
