## Function dependency analysis

Use `npm install` to install dependency for xml parse

Then run `doxygen Doxyfile` to generate analysis xml

Switch to `doc/xml`, run `xsltproc combine.xslt index.xml >all.xml` to merge xml

Switch back to project root, run `node parse.mjs > result.json` to generate dependency json

Run `python3 draw.py > api.json` to generate call graph and API list

After build, run `node parse.2.mjs` to find if any public function is missing

## Prompt and compress

Declare `NEWAPI` or `COMPRESS` in cmake to enable override functions

Any function declared in QuEST_newapi.c or QuEST_compress.c will override function defined in CPU folder, when OVERRIDABLE macro is added to origin functions.

See usage example in .github/workflows/main.yml to build and generate report png.