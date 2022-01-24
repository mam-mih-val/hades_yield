add_executable(yield yield.cc yield.h)
target_link_libraries(yield PUBLIC ROOT::EG at_task_main)

add_executable(super_event super_event.cc super_event.h)
target_link_libraries(super_event PUBLIC ROOT::EG ROOT::Minuit at_task_main)