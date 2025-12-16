#ifndef LOGGER_H
#define LOGGER_H

#include <memory>
#include <mutex>
#include <string>
#include <utility>

#include <spdlog/spdlog.h>

// Thread-safe singleton that configures spdlog for console and rotating file output.
class LogManager {
public:
    static LogManager& Instance();

    void Initialize(const std::string& log_file = "gsc.log",
                    spdlog::level::level_enum level = spdlog::level::info);

    std::shared_ptr<spdlog::logger> Logger();
    void SetLevel(spdlog::level::level_enum level);

    template <typename... Args>
    void Debug(const char* fmt, Args&&... args) {
        Logger()->debug(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Info(const char* fmt, Args&&... args) {
        Logger()->info(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Warn(const char* fmt, Args&&... args) {
        Logger()->warn(fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Error(const char* fmt, Args&&... args) {
        Logger()->error(fmt, std::forward<Args>(args)...);
    }

private:
    LogManager() = default;
    void ConfigureLogger(const std::string& log_file, spdlog::level::level_enum level);

    std::shared_ptr<spdlog::logger> logger_;
    std::once_flag init_flag_;
};

#endif  // LOGGER_H
