#!/bin/bash

# Git快速提交脚本
# 使用方法: ./git-push.sh "提交信息"

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 显示当前状态
echo -e "${YELLOW}=== Git 状态 ===${NC}"
git status --short

# 检查是否有更改
if [ -z "$(git status --porcelain)" ]; then
    echo -e "${GREEN}没有需要提交的更改${NC}"
    exit 0
fi

# 获取提交信息
if [ $# -eq 0 ]; then
    echo -e "${YELLOW}请输入提交信息:${NC}"
    read -r commit_message
else
    commit_message="$1"
fi

# 如果没有输入提交信息，使用默认信息
if [ -z "$commit_message" ]; then
    commit_message="更新代码 $(date '+%Y-%m-%d %H:%M:%S')"
fi

# 显示将要提交的文件
echo -e "\n${YELLOW}=== 将要提交的更改 ===${NC}"
git status --short

# 询问是否继续
echo -e "\n${YELLOW}是否继续提交? (y/n)${NC}"
read -r answer
if [ "$answer" != "y" ] && [ "$answer" != "Y" ]; then
    echo -e "${RED}取消提交${NC}"
    exit 0
fi

# 添加所有更改
echo -e "\n${GREEN}1. 添加所有更改...${NC}"
git add -A

# 提交
echo -e "${GREEN}2. 提交更改...${NC}"
git commit -m "$commit_message"

# 检查提交是否成功
if [ $? -ne 0 ]; then
    echo -e "${RED}提交失败！${NC}"
    exit 1
fi

# 获取当前分支名
current_branch=$(git branch --show-current)
echo -e "${GREEN}当前分支: $current_branch${NC}"

# 推送到远程仓库
echo -e "${GREEN}3. 推送到远程仓库...${NC}"
git push myfork "$current_branch"

# 检查推送是否成功
if [ $? -eq 0 ]; then
    echo -e "\n${GREEN}✓ 成功推送到 myfork/$current_branch！${NC}"
else
    echo -e "\n${RED}推送失败！${NC}"
    echo -e "${YELLOW}尝试使用 --set-upstream 选项...${NC}"
    git push --set-upstream myfork "$current_branch"
    
    if [ $? -eq 0 ]; then
        echo -e "\n${GREEN}✓ 成功推送并设置上游分支！${NC}"
    else
        echo -e "\n${RED}推送仍然失败，请检查网络连接或仓库权限${NC}"
        exit 1
    fi
fi

# 显示最新的提交
echo -e "\n${YELLOW}=== 最新提交 ===${NC}"
git log --oneline -1
