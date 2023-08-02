#ifndef THREAD_POOL_MANAGER_HPP
#define THREAD_POOL_MANAGER_HPP

#include <thread>
#include <future>
#include <queue>
#include <mutex>
#include <tuple>
#include <atomic>
#include <condition_variable>
#include <exception>

struct MultithreadLevelException : public std::exception
{
    const char *what() const throw()
    {
        return "MultithreadLevelException: Cannot start a thread pool with less then 1 thread!";
    }
};

/**
 * @brief The following class defines a thread library whose aim is to create a persistent pool of threads,
 * allowing for the assignement of new tasks without the need to manually create and join each single thread.
 *
 */
template <class ProtectionMethod>
class thread_pool_manager
{
private:
    std::vector<std::thread> threadsPool;
    std::queue<std::function<void()>> tasksQueue;
    std::mutex poolMutex;
    std::mutex queueMutex;
    std::condition_variable poolCondition;
    bool terminatePool;

    /**
     * @brief The following method creates as many threads as specified by the numOfThreads param.
     *
     * @param numOfThreads
     */
    void createThreads(int numOfThreads)
    {
        terminatePool = false;
        if (numOfThreads > 0)
        {
            for (int i = 0; i < numOfThreads; i++)
            {
                threadsPool.push_back(
                    std::thread(&thread_pool_manager::waitForTasks, this));
            }
        }
        else
        {
            throw MultithreadLevelException();
        }
    }

    /**
     * @brief The following method keeps each created thread on idle, waiting for a new task to be executed.
     * When a new task becomes available, one of the idle threads will be awakened in order to process the given task.
     *
     */
    void waitForTasks()
    {
        std::function<void()> taskToPerform;
        while (true)
        {
            {
                std::unique_lock<std::mutex> lock(queueMutex);
                poolCondition.wait(lock, [&]
                                   { return !tasksQueue.empty() || terminatePool; });
                if (!tasksQueue.empty())
                {
                    taskToPerform = tasksQueue.front();
                    tasksQueue.pop();
                }
                else
                {
                    return;
                }
            }
            if (taskToPerform != nullptr)
            {
                taskToPerform();
            }
        }
    }

public:
    ProtectionMethod protection;

    thread_pool_manager(const thread_pool_manager &) = delete;
    thread_pool_manager &operator=(const thread_pool_manager &) = delete;

    /**
     * @brief The following constructor launches the creation of as many threads
     * as the hardware concurrency of the running machine.
     *
     */
    thread_pool_manager()
    {
        createThreads(std::thread::hardware_concurrency());
    }

    /**
     * @brief The following constructor launches the creation of as many threads
     * as specified by the argument.
     *
     * @param numOfThreads
     */
    thread_pool_manager(int numOfThreads)
    {
        createThreads(numOfThreads);
    }

    /**
     * @brief The following method returns the number of created threads.
     *
     * @return int
     */
    int getPoolSize()
    {
        return threadsPool.size();
    }

    /**
     * @brief The following method pushes the new task into the queue of tasks, so to execute it as soon as a thread
     * will be on idle.
     *
     * @param threadIndex
     */
    void executeTask(std::function<void()> &&newTask)
    {
        {
            std::lock_guard<std::mutex> lock(queueMutex);
            tasksQueue.push(std::function<void()>(newTask));
        }
        poolCondition.notify_one();
    }

    /**
     * @brief The following method shutdowns the pool of threads.
     *
     */
    void shutdown()
    {
        {
            std::lock_guard<std::mutex> thMutex(poolMutex);
            terminatePool = true;
        }
        poolCondition.notify_all();

        // Join all threads.
        for (std::thread &actThread : threadsPool)
        {
            actThread.join();
        }
        // Vector and queue emptying
        /* threadsPool.clear();
        std::queue<std::function<void()>>().swap(tasksQueue); */
    }
};

#endif // THREAD_POOL_MANAGER_HPP