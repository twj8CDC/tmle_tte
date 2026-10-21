
data = read_delim(I("
t	t	t	t	t	t	t	t	t	t	t
c	t	t	t	t	t	t	t	t	t	
c	t	t	t	t	t	t	t	t		
c	t	t	t	t	t	t	t	t		
c	c	t	t	t	t	t	t			
c	c	t	t	t	t	t	t			
c	c	t	t	t	t	t	t			
c	c	c	c	t	t	t				
c	c	c	c	c	t					
t	t	t	t							
c	c	c	c	c	c	c	c	c	c	c
c	c	c	c	c	c	c	c			
c	c	c	c	c	c					
c	c	c	c	c						
c	c	c	c							
c	c	c	c							
c	c	c	c							
c	c	c								
c	c	c								
c	c																																											
"
),
delim="\\n",
col_names=FALSE
)

df1 = data %>% 
    rename(x1=X1) %>% 
    mutate(id=row_number()) %>% 
    mutate(tx = str_split(x1, "")) %>%
    unnest(tx) %>% 
    group_by(id)  %>% 
    mutate(t = row_number()-1) %>% 
    mutate(y = ifelse(t==max(t), 1,0)) %>% 
    group_by(id, tx) %>% 
    mutate(dur = row_number())  %>% 
    mutate(dur = ifelse(tx=="c",1, dur)) %>% 
    print(n=20)

df1 %>% tail(10)
lm1 = glm(y ~ tx + I(t^2) + t + dur + I(dur^2), family = binomial, data=bind_rows(df1,df1))    
lm1 = glm(y ~ tx + factor(t) + factor(dur), family = binomial, data=df1)    

summary(lm1)


df2 = tibble(
        id = rep(1:20, each=11), 
        tx1 = c(rep("t", 110), rep("c", 110))
    ) %>% 
    group_by(id) %>% 
    mutate(t = row_number()-1) %>% 
    left_join(df1, by= c("id","t")) %>% 
    mutate(tx2 = ifelse(!is.na(tx), tx, tx1)) %>% 
    group_by(id, tx2) %>% 
    mutate(dur2 = row_number()) %>% 
    ungroup %>% 
    # select(id, tx1, t, dur2) %>% 
    # rename(tx = tx1, dur=dur2) %>%
    select(id, tx2, t, dur2) %>% 
    rename(tx = tx2, dur=dur2) %>%
    mutate(dur = ifelse(tx=="c", 1, dur)) %>% 
    print(n=20)


df2 %>%
    ungroup %>% 
    mutate(preds = predict(lm1, newdata=df2, type="response")) %>% 
    mutate(preds = pmin(0.999, pmax(.0001, preds))) %>% 
    group_by(id) %>% 
    select(id,t,preds) %>% 
    pivot_wider(names_from=t, values_from=preds)


df2 %>%
    ungroup %>% 
    mutate(preds = predict(lm1, newdata=df2, type="response")) %>% 
    mutate(preds = pmin(0.999, pmax(.0001, preds))) %>% 
    group_by(id) %>% 
    mutate(sv = cumprod(1-preds)) %>% 
    ungroup() %>% 
    select(id, t, sv) %>% 
    pivot_wider(names_from =t, values_from=sv) %>% 
    mutate(tx_id = ifelse(id <=10, 1,0)) %>% 
    group_by(tx_id) %>% 
    summarise(across(everything(), mean))


df1 = data %>% 
    rename(x1=X1) %>% 
    mutate(id=row_number()) %>% 
    mutate(x12 = str_pad(x1, side="right", pad="a", width=max(nchar(x1)))) %>%  
    mutate(cens = str_split(x12,"")) %>% 
    unnest(cens) %>% 
    mutate(cens = as.numeric(cens)) %>% 
    # filter(!is.na(spl)) %>% 
    group_by(id) %>% 
    mutate(t = row_number()-3) %>% 
    mutate(tx = ifelse(str_detect(x1, "tx"), 1,0))  %>% 
    filter(t>=0) %>% 
    filter(!is.na(cens))%>%
    select(id, t, tx, cens) %>%
    rename(event=cens) %>% 
    print(n=30)



lm1 = glm(event ~ tx * factor(t), family=binomial, data=df1)

eps=0.00001
df2 = df1 %>% 
    ungroup %>% 
    mutate(
    preds=pmax(predict(lm1, typ="response"),eps)
) %>% 
    group_by(id) %>% 
    mutate(sv = cumprod(1-preds))

df2
df2 %>% filter(t ==2)
df2 %>% filter(preds==max(preds))
df2

cf = bind_rows(df1,df1) %>% 
    ungroup() %>% 
    mutate(preds = 
        c(
            predict(lm1, df1 %>% mutate(tx =1), 'response'),
            predict(lm1, df1 %>% mutate(tx =0), 'response')
        )
    ) %>% 
    mutate(tx = c(rep(1,nrow(df1)), rep(0,nrow(df1)))) %>% 
    group_by(tx, t) %>% 
    summarise(mpreds = mean(preds)) %>% print(n=100) %>% 
    mutate(sv = round(cumprod(1-mpreds),3)) %>% 
    print(n=100)

cf$sv[cf$tx==1] - cf$sv[cf$tx==0]

df2
