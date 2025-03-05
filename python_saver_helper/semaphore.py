#!python
import win32event


def sema_wait(semkeyStr):
    hSemaphore = win32event.OpenSemaphore(win32event.SYNCHRONIZE,False,semkeyStr)
    if hSemaphore == 0:
        print('sema_wait: Unable to open the semaphore handle.')
    else:
        print('sema_wait: waiting')
        wait_failed =  win32event.WaitForSingleObject(hSemaphore,win32event.INFINITE)
        if wait_failed:
          print('sema_wait: Unable to open the semaphore due to failure.')
        else:
          print('sema_wait: wait over')


def sema_post(semkeyStr):
    hSemaphore = win32event.OpenSemaphore(win32event.SEMAPHORE_MODIFY_STATE,false,semkeyStr)
    if hSemaphore == 0:
        print('sema_post: Unable to open the semaphore handle.')
    else:
        hSemaphore = win32event.ReleaseSemaphore(hSemaphore, 1, NULL)
        if hSemaphore == 0:
           print('sema_post: Unable to open the semaphore handle.')
           lastError = win32event.GetLastError();
           errMsg = win32api.FormatMessage(lastError, flags=FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS | FORMAT_MESSAGE_MAX_WIDTH_MASK)
           win32api.CloseHandle(hSemaphore)
           print('sema_post: Unable to post the semaphore with key ' + semkeyStr + ',  error code:' + str(lastError) + ': ' + errMsg)
        else:
           win32api.CloseHandle(hSemaphore);


sema_wait(semkeyStr =  '24537')

