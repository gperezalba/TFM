pragma solidity ^0.4.18;

contract SimpleWallet {

    modifier onlyowner() {
        if (owner == msg.sender){
            _;
        } else {
            revert();
        }
    }

    uint someVar;
    address owner;

    mapping(address => Permission) myAddressMapping;
    
    struct Permission {
        bool isAllowed;
        uint maxAmount;
    }

    event someoneAddedToSenderList(address whoAdded, address whoIsAllowed, uint amountCanSend);

    constructor() public {
        owner = msg.sender;
    }

    function addAddressToSensersList(address permitted, uint maxTransferAmount) public onlyowner {
        myAddressMapping[permitted] = Permission(true, maxTransferAmount);
        emit someoneAddedToSenderList(msg.sender, permitted, maxTransferAmount);
    }

    function removeAddressToSendersList(address forbidded) public onlyowner {
        delete myAddressMapping[forbidded]; //or myAddressMapping[forbidded] = Permission(false, 0);
    }

    function sendFunds(address receiver, uint amountInWei) public {
        if(myAddressMapping[msg.sender].isAllowed){
            if(myAddressMapping[msg.sender].maxAmount >= amountInWei){
                receiver.transfer(amountInWei);
            } else {
                revert();
            }
        } else {
            revert();
        }
    }
    
    function () public payable {

    }
}