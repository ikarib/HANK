CC = nvcc
NAME = estimate
DFLS = DFLS
CFLAGS = -O2 -I${DFLS}
LDFLAGS = -L${DFLS}

default: ${NAME}

clean:
	rm -f ${NAME}

${NAME}: ${NAME}.cu
	${CC} $(CFLAGS) ${LDFLAGS} $< -o $@ -lnewuoa_h

test: ${NAME}.out ${NAME}-K80.out
	@if diff $^ >${NAME}-test.diff; then \
	    msg="test of C version successful"; \
	    /bin/echo -e "\e[1;32m$${msg}\e[0;39;49m"; \
	else \
	    msg="test of C version failed (see file ${NAME}-test.diff)"; \
	    /bin/echo -e "\e[1;31m$${msg}\e[0;39;49m"; \
	fi
