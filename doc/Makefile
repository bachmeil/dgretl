all: index base

index:
	pandoc -s index.md -o index.html --template=template.html

base:
	pandoc -s doc-base.md -o doc-base.html --template=template.html
