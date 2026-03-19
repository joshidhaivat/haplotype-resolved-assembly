#ifndef THREAD_SAFE_QUEUE_H
#define THREAD_SAFE_QUEUE_H

#include <queue>
#include <mutex>
#include <condition_variable>
#include <chrono>

template<typename T>
class ThreadSafeQueue {
private:
    std::queue<T> queue_;
    mutable std::mutex mutex_;
    std::condition_variable cond_not_empty_;
    std::condition_variable cond_not_full_;
    bool finished_ = false;
    size_t max_size_;  // Maximum queue size for backpressure

public:
    // Constructor with optional max size (0 = unlimited)
    explicit ThreadSafeQueue(size_t max_size = 0) : max_size_(max_size) {}

    void push(T value) {
        std::unique_lock<std::mutex> lock(mutex_);
        
        // If max_size is set, wait if queue is full (backpressure)
        if (max_size_ > 0) {
            cond_not_full_.wait(lock, [this] { 
                return queue_.size() < max_size_ || finished_; 
            });
        }
        
        queue_.push(std::move(value));
        lock.unlock();
        cond_not_empty_.notify_one();
    }

    bool pop(T& value) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_not_empty_.wait(lock, [this] { return !queue_.empty() || finished_; });
        
        if (queue_.empty()) return false;
        
        value = std::move(queue_.front());
        queue_.pop();
        
        lock.unlock();
        cond_not_full_.notify_one();  // Notify producer that space is available
        return true;
    }

    // Non-blocking try_pop with timeout
    bool try_pop(T& value, std::chrono::milliseconds timeout = std::chrono::milliseconds(100)) {
        std::unique_lock<std::mutex> lock(mutex_);
        
        if (!cond_not_empty_.wait_for(lock, timeout, [this] { return !queue_.empty() || finished_; })) {
            return false;  // Timeout
        }
        
        if (queue_.empty()) return false;
        
        value = std::move(queue_.front());
        queue_.pop();
        
        lock.unlock();
        cond_not_full_.notify_one();
        return true;
    }

    void finish() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            finished_ = true;
        }
        cond_not_empty_.notify_all();
        cond_not_full_.notify_all();
    }

    size_t size() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.size();
    }

    bool is_finished() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return finished_;
    }
    
    bool empty() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }
};

#endif // THREAD_SAFE_QUEUE_H