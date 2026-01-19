#pragma once

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <type_traits>
#include <vector>

namespace eastr {

/**
 * Simple thread pool for parallel task execution.
 */
class ThreadPool {
public:
    /**
     * Create a thread pool with the specified number of threads.
     * @param num_threads Number of worker threads (0 = hardware concurrency)
     */
    explicit ThreadPool(size_t num_threads = 0);

    ~ThreadPool();

    // Non-copyable
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /**
     * Submit a task to the thread pool.
     * @param f Function to execute
     * @param args Arguments to pass to the function
     * @return Future for the result
     */
    template<typename F, typename... Args>
    auto submit(F&& f, Args&&... args)
        -> std::future<std::invoke_result_t<F, Args...>> {

        using return_type = std::invoke_result_t<F, Args...>;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> result = task->get_future();

        {
            std::unique_lock<std::mutex> lock(mutex_);
            tasks_.emplace([task]() { (*task)(); });
            pending_tasks_++;
        }

        condition_.notify_one();
        return result;
    }

    /**
     * Wait for all submitted tasks to complete.
     */
    void wait_all();

    /**
     * Get number of worker threads.
     */
    size_t num_threads() const { return workers_.size(); }

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;

    std::mutex mutex_;
    std::condition_variable condition_;
    std::condition_variable done_condition_;

    std::atomic<size_t> active_tasks_{0};
    std::atomic<size_t> pending_tasks_{0};
    bool stop_{false};
};

} // namespace eastr
