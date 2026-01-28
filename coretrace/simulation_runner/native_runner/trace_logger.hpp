#ifndef SOLTRACE_TRACE_LOGGER_H
#define SOLTRACE_TRACE_LOGGER_H

#include <memory>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

namespace SolTrace::NativeRunner
{

    class TraceLogger
    {
    public:
        TraceLogger();
        ~TraceLogger();

        void error_log(const std::string &msg);
        void print_log(std::ostream &os) const;

    private:
        mutable std::mutex message_mutex;
        std::vector<std::string> messages;
    };

    using trace_logger_ptr = std::shared_ptr<TraceLogger>;

    template <typename... Args>
    inline auto make_trace_logger(Args &&...args)
    {
        return std::make_shared<TraceLogger>(std::forward<Args>(args)...);
    }

} // namespace SolTrace::NativeRunner

#endif
