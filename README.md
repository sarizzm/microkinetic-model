# microkinetic-model
microkinetic model of chemical reaction,  微观动力学模型总结

书籍进行更新，是我以前看的最新版，也可以到作者管网去看看有没有再次更新

https://space.bilibili.com/80446042/channel/seriesdetail?sid=2463662&ctype=0

# 1 介绍：

其他信息看完版，最好把catmap 官方说明文档自己看一遍，尤其是自己相关的。

```
https://github.com/SUNCAT-Center/catmap
```

这里面也有一些书籍和示例供大家阅读：

```
https://github.com/sarizzm/microkinetic-model
https://space.bilibili.com/80446042/channel/seriesdetail?sid=2463662&ctype=0
```

# 2 主要作用：

**介绍一些奇技淫巧：**

## 2.1 结构书写

1. **分子式换一种写法**

```
# H2O 可以写成 H2O,OH2,HOH,等，不影响catmap 计算，但是不要加奇怪的前缀和后缀
# 比如：line-H2O, 因为解析的过程 catmap 无法识别line是啥东西，报错
```

1. 同一种分子在不同位点，可以分子式一样，会自动识别为不同分子
2. 如果为气态就不要按照第一种方法，因为气态会调用ase进行识别，ase识别不出来分子也会报错，同时可以借助catmap提供接口，进行相关气态的补充。

2.2 catmap 源码讲解（有时间的话）





欢迎加微信`xgdzzm`，因为微信群只有七天有效期，拉你入群讨论学习，愿你结交更多学友，不会尽量不要问我，因为我会的不多。

![](https://github.com/sarizzm/microkinetic-model/blob/main/picture/wx.jpg?raw=true)

