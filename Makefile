# Makefile
# Makefile for package checking building, installing, uninstalling, etc.

PKG=ebimetagenomics
VERSION=0.2

DEFAULT:
	make check
	make install

check:
	R CMD check pkg

check-cran:
	R CMD check --as-cran pkg

build:
	R CMD build pkg

install:
	make build
	R CMD INSTALL $(PKG)_$(VERSION).tar.gz

install-pkg:
	R CMD INSTALL pkg

remove:
	R CMD REMOVE $(PKG)

clean:
	rm -rf *~ $(PKG)_*.tar.gz pkg.Rcheck


update:
	svn update
	svn log|less

commit:
	svn commit
	make update



# eof

