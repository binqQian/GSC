#include "logger.h"

#include <iostream>

#include <spdlog/sinks/dist_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace {
constexpr std::size_t kMaxLogFileSize = 10 * 1024 * 1024;  // 10MB
constexpr std::size_t kMaxRotatedFiles = 3;
}  // namespace

LogManager& LogManager::Instance() {
    static LogManager instance;
    return instance;
}

void LogManager::Initialize(const std::string& log_file, spdlog::level::level_enum level) {
    std::call_once(init_flag_, [this, log_file, level]() {
        try {
            ConfigureLogger(log_file, level);
        } catch (const spdlog::spdlog_ex& ex) {
            std::cerr << "Failed to initialize logger: " << ex.what() << std::endl;
        }
    });
}

std::shared_ptr<spdlog::logger> LogManager::Logger() {
    Initialize();
    if (logger_) {
        return logger_;
    }
    return spdlog::default_logger();
}

void LogManager::SetLevel(spdlog::level::level_enum level) {
    auto log = Logger();
    if (log) {
        log->set_level(level);
        spdlog::set_level(level);
    }
}

void LogManager::ConfigureLogger(const std::string& log_file, spdlog::level::level_enum level) {
    auto dist_sink = std::make_shared<spdlog::sinks::dist_sink_mt>();

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_pattern("[%H:%M:%S.%e] [%^%l%$] %v");

    auto rotating_sink =
        std::make_shared<spdlog::sinks::rotating_file_sink_mt>(log_file, kMaxLogFileSize, kMaxRotatedFiles);
    rotating_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");

    dist_sink->add_sink(console_sink);
    dist_sink->add_sink(rotating_sink);

    logger_ = std::make_shared<spdlog::logger>("gsc_logger", dist_sink);
    logger_->set_level(level);
    logger_->flush_on(spdlog::level::info);

    spdlog::set_default_logger(logger_);
    spdlog::set_level(level);
}
