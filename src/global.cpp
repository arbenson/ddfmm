// TODO (Austin): This needs to be compatible with the Elemental license

#ifndef RELEASE
// If we are not in RELEASE mode, then implement wrappers for a CallStack
#include <stack>
#include <iostream>
#include <sstream>

std::stack<std::string> callStack;

void PushCallStack( std::string s )
{ 
#ifdef HAVE_OPENMP
    if( omp_get_thread_num() != 0 )
        return;
#endif // HAVE_OPENMP
    ::callStack.push(s); 
}

void PopCallStack()
{ 
#ifdef HAVE_OPENMP
    if( omp_get_thread_num() != 0 )
        return;
#endif // HAVE_OPENMP
    ::callStack.pop(); 
}

void DumpCallStack( std::ostream& os )
{
    std::ostringstream msg;
    while( ! ::callStack.empty() ) {
        msg << "[" << ::callStack.size() << "]: " << ::callStack.top() << "\n";
        ::callStack.pop();
    }
    os << msg.str();
    os.flush();
}
#endif  // RELEASE
