include makefile.${HOST}

LIB_SRC	=	wave3d.cpp	kernel3d.cpp	mlib3d.cpp	wave3d_setup.cpp	wave3d_eval.cpp	wave3d_check.cpp \
		vecmatop.cpp	parallel.cpp

LIB_OBJ	 = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	tt.cpp ss.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libwave.a:	${LIB_OBJ}
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

tt:	libwave.a tt.o
	${CXX} -o tt tt.o libwave.a ${LDFLAGS}

ss: ss.o
	${CXX} -o ss ss.o ${LDFLAGS}

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.a

tags:
	etags *.hpp *.cpp

