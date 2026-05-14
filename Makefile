PKG     := samjr
PKGDIR  := $(CURDIR)/$(PKG)
VERSION := $(shell awk '/^Version:/ {print $$2}' $(PKGDIR)/DESCRIPTION)
TARBALL := $(PKG)_$(VERSION).tar.gz

R       ?= R
ROXY    := cd $(PKGDIR) && $(R) --quiet --no-save -e 'roxygen2::roxygenize(".")'

R_FILES    := $(wildcard $(PKGDIR)/R/*.R)
RD_FILES   := $(wildcard $(PKGDIR)/man/*.Rd)
DATA_FILES := $(wildcard $(PKGDIR)/data/*.rda)

TESTMORE_DIR     := $(CURDIR)/testmore
TESTMORE_DIRS    := $(patsubst $(TESTMORE_DIR)/%/script.R,%,$(wildcard $(TESTMORE_DIR)/*/script.R))
TESTMORE_TIMEOUT ?= 300
TESTMORE_STATUS  := $(TESTMORE_DIR)/.status

.PHONY: all doc install build check clean data help testmore testmore-prep testmore-summary $(addprefix testmore-,$(TESTMORE_DIRS))

all: install

help:
	@echo "Targets:"
	@echo "  doc       regenerate $(PKG)/man/*.Rd and NAMESPACE from roxygen comments"
	@echo "  install   doc + R CMD INSTALL $(PKG)"
	@echo "  build     doc + R CMD build (produces $(TARBALL))"
	@echo "  check     build + R CMD check --as-cran on the tarball"
	@echo "  data      rebuild $(PKG)/data/nscod*.rda from testmore/nscod"
	@echo "  testmore  run every testmore/<dir>/script.R and print OK / FAIL"
	@echo "            (use 'make -j N testmore' to run in parallel)"
	@echo "  clean     remove generated tarballs, check dirs, testmore artefacts"

doc: $(R_FILES)
	$(ROXY)

install: doc
	$(R) CMD INSTALL $(PKGDIR)

build: doc
	$(R) CMD build $(PKGDIR)

check: build
	_R_CHECK_SYSTEM_CLOCK_=FALSE $(R) CMD check --as-cran $(TARBALL)

data: install
	cd $(PKGDIR) && $(R) --quiet --no-save -f tools/build-nscod-data.R

# ---------------------------------------------------------------------------
# testmore: run every testmore/<dir>/script.R, print OK / FAIL per test, then
# a summary. Per-test results are recorded in $(TESTMORE_STATUS)/<dir> so the
# summary works under `make -j`.

testmore-prep:
	@rm -rf $(TESTMORE_STATUS)
	@mkdir -p $(TESTMORE_STATUS)

testmore: testmore-prep $(addprefix testmore-,$(TESTMORE_DIRS)) testmore-summary

$(addprefix testmore-,$(TESTMORE_DIRS)): testmore-%: testmore-prep
	@d=$(TESTMORE_DIR)/$*; \
	 cd $$d && rm -f res.out errlog.txt; \
	 t0=$$(date +%s.%N); \
	 timeout $(TESTMORE_TIMEOUT) $(R) --slave --vanilla -e 'source("script.R")' \
	   > /dev/null 2> errlog.txt; rc=$$?; \
	 t1=$$(date +%s.%N); \
	 el=$$(awk -v a=$$t0 -v b=$$t1 'BEGIN{printf "%.1f", b-a}'); \
	 if [ $$rc -ne 0 ] || [ ! -f res.out ]; then \
	   printf "  %-25s %5ss  FAIL  (run error)\n" "$*" "$$el"; \
	   echo "FAIL $$el" > $(TESTMORE_STATUS)/$*; \
	 elif diff --strip-trailing-cr -q res.out res.EXP > /dev/null 2>&1; then \
	   printf "  %-25s %5ss  OK\n" "$*" "$$el"; \
	   echo "OK $$el" > $(TESTMORE_STATUS)/$*; \
	 else \
	   printf "  %-25s %5ss  FAIL  (diff)\n" "$*" "$$el"; \
	   echo "FAIL $$el" > $(TESTMORE_STATUS)/$*; \
	 fi

testmore-summary: $(addprefix testmore-,$(TESTMORE_DIRS))
	@ok=$$(awk '$$1=="OK"'   $(TESTMORE_STATUS)/* 2>/dev/null | wc -l); \
	 fail=$$(awk '$$1=="FAIL"' $(TESTMORE_STATUS)/* 2>/dev/null | wc -l); \
	 tot=$$((ok + fail)); \
	 totT=$$(awk '{s+=$$2} END{printf "%.1f", s}' $(TESTMORE_STATUS)/* 2>/dev/null); \
	 echo "===== testmore: $$ok OK, $$fail FAIL (of $$tot) - $${totT}s ====="; \
	 rm -rf $(TESTMORE_STATUS); \
	 [ $$fail -eq 0 ]

clean:
	rm -f $(CURDIR)/$(TARBALL)
	rm -rf $(CURDIR)/$(PKG).Rcheck
	rm -rf $(TESTMORE_STATUS)
	find $(TESTMORE_DIR) -mindepth 2 -maxdepth 2 -name res.out  -delete
	find $(TESTMORE_DIR) -mindepth 2 -maxdepth 2 -name errlog.txt -delete
