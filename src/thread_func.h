#include <iostream>

template <bool perturb, formula f, color c>
void mariani_silver(const double&, const double&, const mariani_data, const perturb_data* const);

template <bool perturb, formula f, color c>
void initial_perimeter(const double&, const double&, const mariani_data&, const perturb_data* const);


struct mariani_cfg{
        double dy;
        double dx;
        mariani_data m_d;
        perturb_data* p_d;
};

mariani_cfg cfg_pool[25];
node *node_pool[25];

void threadfunc(taskqueue *q){
        while (1){
                node* n = q->dequeue();
                if (n == NULL) { break;}
                n->func(n->ptr);
        }
}

template <bool perturb, formula f, color c>
static void mariani_wrapper(void *ptr){
        mariani_cfg* m_cfg = (mariani_cfg*)ptr;
        initial_perimeter<perturb, f, c>(m_cfg->dy, m_cfg->dx, m_cfg->m_d, m_cfg->p_d);
        mariani_silver<perturb, f, c>(m_cfg->dy, m_cfg->dx, m_cfg->m_d, m_cfg->p_d);
        q.mtx.lock();
        q.current_tasks--;
        if (q.current_tasks == 0){
                q.cv.notify_all();
        }
        q.mtx.unlock();
}
