#include "trace_logger.hpp"

#include <ostream>

namespace SolTrace::NativeRunner
{
    TraceLogger::TraceLogger()
    {
        return;
    }

    TraceLogger::~TraceLogger()
    {
        this->messages.clear();
        return;
    }

    void TraceLogger::error_log(const std::string &msg)
    {
        std::lock_guard<std::mutex> lk(this->message_mutex);
        this->messages.push_back(msg);
        return;
    }
    
    void TraceLogger::print_log(std::ostream &os) const
    {
        std::lock_guard<std::mutex> lk(this->message_mutex);
        for (auto iter = this->messages.cbegin();
             iter != this->messages.cend();
             ++iter)
        {
            os << *iter;
        }
        return;
    }

} // namespace SolTrace::NativeRunner
