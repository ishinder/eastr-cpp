#include "eastr/thread_pool.hpp"

namespace eastr {

ThreadPool::ThreadPool(size_t num_threads) {
    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 1;
    }

    workers_.reserve(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
        workers_.emplace_back([this] {
            while (true) {
                std::function<void()> task;

                {
                    std::unique_lock<std::mutex> lock(mutex_);
                    condition_.wait(lock, [this] {
                        return stop_ || !tasks_.empty();
                    });

                    if (stop_ && tasks_.empty()) {
                        return;
                    }

                    task = std::move(tasks_.front());
                    tasks_.pop();
                    active_tasks_++;
                }

                task();

                {
                    std::unique_lock<std::mutex> lock(mutex_);
                    active_tasks_--;
                    pending_tasks_--;
                    if (pending_tasks_ == 0 && active_tasks_ == 0) {
                        done_condition_.notify_all();
                    }
                }
            }
        });
    }
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(mutex_);
        stop_ = true;
    }
    condition_.notify_all();

    for (auto& worker : workers_) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

void ThreadPool::wait_all() {
    std::unique_lock<std::mutex> lock(mutex_);
    done_condition_.wait(lock, [this] {
        return pending_tasks_ == 0 && active_tasks_ == 0;
    });
}

} // namespace eastr
