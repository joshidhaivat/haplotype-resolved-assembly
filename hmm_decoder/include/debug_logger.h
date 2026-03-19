#ifndef DEBUG_LOGGER_H
#define DEBUG_LOGGER_H

#include <string>
#include <fstream>
#include <mutex>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <thread>

// ============================================================================
// COMPILE-TIME TOGGLE
// ============================================================================
// Comment out or set to 0 to completely disable logging (zero runtime cost)
#define DEBUG_LOGGING_ENABLED 1

// ============================================================================
// LOGGING MACROS
// ============================================================================
#if DEBUG_LOGGING_ENABLED

    // Basic log with automatic file/line/function info
    #define LOG_DEBUG(msg) \
        DebugLogger::instance().log(__FILE__, __LINE__, __func__, msg)

    // Log with stream syntax: LOG_STREAM("value=" << x << " count=" << count)
    #define LOG_STREAM(expr) \
        do { \
            std::ostringstream _log_ss; \
            _log_ss << expr; \
            DebugLogger::instance().log(__FILE__, __LINE__, __func__, _log_ss.str()); \
        } while(0)

    // Log variable name and value: LOG_VAR(myVariable) -> "myVariable = 42"
    #define LOG_VAR(var) \
        do { \
            std::ostringstream _log_ss; \
            _log_ss << #var << " = " << (var); \
            DebugLogger::instance().log(__FILE__, __LINE__, __func__, _log_ss.str()); \
        } while(0)

    // Log vector contents
    #define LOG_VECTOR(vec) \
        do { \
            std::ostringstream _log_ss; \
            _log_ss << #vec << " [size=" << (vec).size() << "]: "; \
            for (size_t _i = 0; _i < std::min((vec).size(), (size_t)10); ++_i) { \
                if (_i > 0) _log_ss << ", "; \
                _log_ss << (vec)[_i]; \
            } \
            if ((vec).size() > 10) _log_ss << " ..."; \
            DebugLogger::instance().log(__FILE__, __LINE__, __func__, _log_ss.str()); \
        } while(0)

    // Conditional logging
    #define LOG_IF(condition, msg) \
        do { if (condition) LOG_DEBUG(msg); } while(0)

    // Function entry/exit tracing
    #define LOG_FUNC_ENTER() LOG_DEBUG(">>> ENTER")
    #define LOG_FUNC_EXIT()  LOG_DEBUG("<<< EXIT")

    // Initialize logger (call once in main)
    #define LOG_INIT(filename) DebugLogger::instance().init(filename)
    
    // Shutdown logger (optional, call at end of main)
    #define LOG_SHUTDOWN() DebugLogger::instance().shutdown()

#else
    // When disabled, all macros expand to nothing (zero cost)
    #define LOG_DEBUG(msg)         ((void)0)
    #define LOG_STREAM(expr)       ((void)0)
    #define LOG_VAR(var)           ((void)0)
    #define LOG_VECTOR(vec)        ((void)0)
    #define LOG_IF(condition, msg) ((void)0)
    #define LOG_FUNC_ENTER()       ((void)0)
    #define LOG_FUNC_EXIT()        ((void)0)
    #define LOG_INIT(filename)     ((void)0)
    #define LOG_SHUTDOWN()         ((void)0)
#endif

// ============================================================================
// LOGGER CLASS (Singleton, Thread-Safe)
// ============================================================================
class DebugLogger {
public:
    // Get singleton instance
    static DebugLogger& instance() {
        static DebugLogger logger;
        return logger;
    }

    // Initialize with log file path (call once at startup)
    void init(const std::string& filename) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (file_.is_open()) {
            file_.close();
        }
        file_.open(filename, std::ios::out | std::ios::trunc);
        if (!file_.is_open()) {
            throw std::runtime_error("Cannot open debug log file: " + filename);
        }
        enabled_ = true;
        start_time_ = std::chrono::steady_clock::now();
        
        // Write header
        file_ << "=== Debug Log Started ===" << std::endl;
        file_ << "Timestamp (ms) | Thread | File:Line | Function | Message" << std::endl;
        file_ << std::string(80, '-') << std::endl;
        file_.flush();
    }

    // Log a message (thread-safe)
    void log(const char* file, int line, const char* func, const std::string& msg) {
        if (!enabled_) return;
        
        std::lock_guard<std::mutex> lock(mutex_);
        if (!file_.is_open()) return;

        // Calculate elapsed time in milliseconds
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_).count();

        // Get thread ID (shortened for readability)
        std::ostringstream tid_ss;
        tid_ss << std::this_thread::get_id();
        std::string tid = tid_ss.str();
        if (tid.length() > 6) tid = tid.substr(tid.length() - 6);

        // Extract just the filename (not full path)
        std::string filename(file);
        size_t pos = filename.find_last_of("/\\");
        if (pos != std::string::npos) {
            filename = filename.substr(pos + 1);
        }

        // Write log entry
        file_ << std::setw(10) << elapsed << " | "
              << std::setw(6) << tid << " | "
              << std::setw(20) << (filename + ":" + std::to_string(line)) << " | "
              << std::setw(20) << func << " | "
              << msg << std::endl;
        
        // Flush immediately for debugging (can be disabled for performance)
        file_.flush();
    }

    // Runtime enable/disable (for conditional logging sections)
    void set_enabled(bool enabled) {
        std::lock_guard<std::mutex> lock(mutex_);
        enabled_ = enabled;
    }

    bool is_enabled() const {
        return enabled_;
    }

    // Clean shutdown
    void shutdown() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (file_.is_open()) {
            file_ << std::string(80, '-') << std::endl;
            file_ << "=== Debug Log Ended ===" << std::endl;
            file_.close();
        }
        enabled_ = false;
    }

    // Destructor ensures file is closed
    ~DebugLogger() {
        shutdown();
    }

private:
    DebugLogger() : enabled_(false) {}
    DebugLogger(const DebugLogger&) = delete;
    DebugLogger& operator=(const DebugLogger&) = delete;

    std::ofstream file_;
    std::mutex mutex_;
    bool enabled_;
    std::chrono::steady_clock::time_point start_time_;
};

#endif // DEBUG_LOGGER_H
