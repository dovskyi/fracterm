#ifndef _QUEUE
#define _QUEUE
#include <iostream>
#include <mutex>
#include <condition_variable>

struct node{
        void (*func)(void*);
        void *ptr;
        node *next;
};

struct taskqueue{
        node* head = NULL;
        node* tail = NULL;

        std::mutex mtx;
        std::condition_variable cv;
        int current_tasks = 0;
        bool stop = 0;

        void enqueue(node* n){
                //allocation is handled outside
                mtx.lock();
                current_tasks++;
                n->next = NULL;
                if (head == NULL){
                        head = n;
                        tail = n;
                        cv.notify_one();
                        mtx.unlock();
                        return;
                }
                tail->next = n;
                tail = n;
                cv.notify_one();
                mtx.unlock();
                return;
        }

        node* dequeue(){
                std::unique_lock<std::mutex> lock(mtx);
                
                while (head == NULL && !stop){
                        cv.wait(lock);
                }
                
                if (head == NULL) { return NULL;}

                node* temp = head;
                head = head->next;
                //deallocation is handled outside
                return temp;
        }

        void wait4threads(){
                std::unique_lock<std::mutex> lock(mtx);
                while (current_tasks != 0) {
                        cv.wait(lock); 
                }
        }
};

#endif
