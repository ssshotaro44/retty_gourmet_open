with open('urls.txt') as f:
    for url in f:
        url_list = url.replace('\n', '').split('/')
        print(url_list)