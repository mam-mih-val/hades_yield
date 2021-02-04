add_executable(yield yield.cc yield.h)
target_link_libraries(yield PUBLIC ROOT::EG at_task_main)